#include <stddef.h>

#include "bwt_internal.h"
#include "unittests.h"

struct cstr_li_durbin_preproc
{
    cstr_alphabet alpha;
    cstr_suffix_array *sa;
    struct c_table *ctab;
    struct o_table *otab;
    struct o_table *rotab;
};

void cstr_free_li_durbin_preproc(cstr_li_durbin_preproc *preproc)
{
    free(preproc->sa);
    free(preproc->ctab);
    free(preproc->otab);
    free(preproc->rotab);
    free(preproc);
}

cstr_li_durbin_preproc *
cstr_li_durbin_preprocess(cstr_const_sslice x)
{
    cstr_li_durbin_preproc *preproc = cstr_malloc(sizeof *preproc);

    cstr_init_alphabet(&preproc->alpha, x);

    // Map the string into an integer slice so we can build the suffix array.
    cstr_uislice *u_buf = cstr_alloc_uislice(x.len);
    cstr_alphabet_map_to_uint(*u_buf, x, &preproc->alpha);
    cstr_const_uislice u = CSTR_SLICE_CONST_CAST(*u_buf);

    cstr_sslice *w_buf = cstr_alloc_sslice(x.len); // Mapping x into w
    cstr_alphabet_map(*w_buf, x, &preproc->alpha);
    cstr_const_sslice w = CSTR_SLICE_CONST_CAST(*w_buf);

    preproc->sa = cstr_alloc_uislice(x.len); // Get the suffix array

    cstr_sslice *bwt_buf = cstr_alloc_sslice(x.len); // Get a buffer for the the bwt string
    cstr_const_sslice bwt = CSTR_SLICE_CONST_CAST(*bwt_buf);

    // Then build the suffix arrays. Build the reverse first,
    // so we don't need a copy for it later. We need to remember
    // the forward one, but not the backward one.
    // The PREFIX stuff is so we only reverse the prefix up to the sentinel,
    // but we do not include the sentinel in the reversal.
    CSTR_REV_SLICE(CSTR_PREFIX(*u_buf, -1)); // Both of these are representations of the input string.
    CSTR_REV_SLICE(CSTR_PREFIX(*w_buf, -1)); // Just with different types. We reverse to build RO

    cstr_sais(*preproc->sa, u, &preproc->alpha);
    cstr_bwt(*bwt_buf, w, *preproc->sa);
    preproc->ctab = cstr_build_c_table(bwt, preproc->alpha.size);
    preproc->rotab = cstr_build_o_table(bwt, preproc->ctab);

    // The C table is the same in either direction, but we need
    // to rebuild the suffix array and BWT to get the O table
    CSTR_REV_SLICE(CSTR_PREFIX(*u_buf, -1)); // Reverse the input back to the forward direction
    CSTR_REV_SLICE(CSTR_PREFIX(*w_buf, -1)); // so we can build the forward BWT and O table

    cstr_sais(*preproc->sa, u, &preproc->alpha);
    cstr_bwt(*bwt_buf, w, *preproc->sa);
    preproc->otab = cstr_build_o_table(bwt, preproc->ctab);

    // We don't need the mapped string nor the BWT any more.
    // The information we need is all in the tables in preproc.
    free(bwt_buf);
    free(w_buf);
    free(u_buf);

    return preproc;
}

// MARK: Cigar stuff

static inline long long scan_next(long long i, const char *edits)
{
    long long j = i;
    for (; j >= 0; j--)
    {
        if (edits[i] != edits[j])
        {
            return j;
        }
    }
    // If we get here, we scanned the last segment and they are all equal,
    // so we should return -1 as past the last (in reverse order)
    return -1;
}

static void edits_to_cigar(char *cigar_buf, const char *edits)
{
    long long edits_len = (long long)strlen(edits);
    for (long long i = edits_len - 1, j = scan_next(i, edits); i >= 0; i = j, j = scan_next(i, edits))
    {
        cigar_buf += sprintf(cigar_buf, "%d%c", (int)(i - j), (char)edits[i]);
    }
    *cigar_buf = '\0';
}

// MARK: State machine for searching
enum states
{
    EMIT_,
    DONE_,
    REC_,
    MATCH_,
    INSERT_,
    DELETE_
};

struct state
{
    enum states state;
    union
    {
        long long left; // LEFT for search
        long long next; // NEXT for emit
    };
    union
    {
        long long right; // RIGHT for search
        long long end;   // END for emit
    };
    long long pos; // Index into p
    long long d;   // Number of edits left
    uint8_t a;     // Current letter we are "deleting"
    char *edits;   // Edit operations so far, when searching
};

// macros for creating new states
#define DONE() ((struct state){.state = DONE_})
#define EMIT(NEXT, END, CIGAR) \
    ((struct state){.state = EMIT_, .next = NEXT, .end = END, .cigar = CIGAR})
#define REC(LEFT, RIGHT, POS, D, EDITS) \
    ((struct state){.state = REC_, .left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS})
#define MATCH(LEFT, RIGHT, POS, D, EDITS, A) \
    ((struct state){.state = MATCH_, .left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS, .a = A})
#define INSERT(LEFT, RIGHT, POS, D, EDITS) \
    ((struct state){.state = INSERT_, .left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS})
#define DELETE(LEFT, RIGHT, POS, D, EDITS, A) \
    ((struct state){.state = DELETE_, .left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS, .a = A})

// MARK: Stack used for CPS.
// Quite simple implementation, but it suffices for what we need here.

// FIXME: I should be able to work out the max stack depth a priori so I don't have to
// waste time on resizing...
struct stack
{
    size_t size;
    size_t used;
    struct state frames[];
};

// We only shrink a stack when it is 1/4 used, and then only to 1/2 size,
// so memory we pop off is still available until the next stack action.
#define MIN_STACK_SIZE 256
static struct stack *new_stack(void)
{
    struct stack *stack = cstr_malloc_header_array(offsetof(struct stack, frames),
                                                   sizeof stack->frames[0], MIN_STACK_SIZE);
    stack->size = MIN_STACK_SIZE;
    stack->used = 0;
    return stack;
}

static inline struct stack *resize_stack(struct stack *stack)
{
    return stack;
}

static inline void push_frame(struct stack **stack, struct state s)
{
    if ((*stack)->used == (*stack)->size)
    {
        (*stack)->size *= 2;
        *stack = cstr_realloc_header_array(stack,
                                         offsetof(struct stack, frames),
                                         sizeof (*stack)->frames[0],
                                         (*stack)->size);
    }
    (*stack)->frames[(*stack)->used++] = s;
}

static inline struct state pop_frame(struct stack **stack)
{
    return (*stack)->frames[--(*stack)->used];
}

// MARK: iterator for searches
// If I implement other approximative matcheres, then
// this should be the Li & Durbin specific, and we need a
// vtab as for the exact matcher to distinguish between
// them.
struct cstr_approx_matcher
{
    struct state current_state;
    cstr_sslice *p_buf;
    cstr_const_sslice p;
    cstr_li_durbin_preproc *preproc;
    // cstr_sslice_buf *edits_buf;
    char *edits_buf;
    char *cigar_buf;
    struct stack *stack;
};

#define S itr->current_state
#define PUSH(STATE) push_frame(&itr->stack, STATE)
#define POP() pop_frame(&itr->stack)

cstr_approx_match cstr_approx_next_match(cstr_approx_matcher *itr)
{
    // Get these into our name space so we can use the macros
    struct c_table *ctab = itr->preproc->ctab;
    struct o_table *otab = itr->preproc->otab;

    // When called, dispatch on the current state...
    goto dispatch;

pop_next: // === Pop the next state from the stack and dispatch it =======================================
    S = POP();
    // fall through...

dispatch: // === Jump to the state in current_state ======================================================
    // clang-format off
    switch (S.state) {
        case EMIT_:   goto emit;
        case DONE_:   goto done;
        case REC_:    goto rec;
        case MATCH_:  goto match;
        case INSERT_: goto insert;
        case DELETE_: goto delete;
    }
    // clang-format on

rec: // === Processing initial recursion states ==========================================================
    if (S.left == S.right || S.d < 0)
    {
        // No matches here, so continue with the next continuation...
        goto pop_next;
    }

    if (S.pos < 0)
    {
        // We have a match!
        *S.edits = '\0'; // terminate this string of edits
        edits_to_cigar(itr->cigar_buf, itr->edits_buf);

        // ---- Go for emit ------------------------
        S.next = S.left;
        S.end = S.right;
        goto emit;
    }

    // Otherwise, continue the recursion with the next match operation.
    // The match operation starts with a == 1 because we don't want
    // the sentinel.
    // ---- Go for match/subs (with the first letter) -----------
    S.a = 1;
    goto match;

match: // === Processing matching states ==============================================================
    if (S.a == itr->preproc->alpha.size)
    {
        // We are through the alphabet, so there is nothing more to match.
        // Try inserting instead.
        goto insert;
    }

    // === Recurse, then continue with the next match afterwards =====================================
    // --- First we put the match continuation on the stack...   -------------------------------------
    PUSH(MATCH(S.left, S.right, S.pos, S.d, S.edits, (uint8_t)(S.a + 1)));

    // --- Then we continue with the recursive operation on pos-1 ------------------------------------
    S.left = C(S.a) + O(S.a, S.left);
    S.right = C(S.a) + O(S.a, S.right);
    S.d -= itr->p.buf[S.pos] != S.a;
    S.pos--;
    *S.edits++ = 'M';
    goto rec;

insert: // === Processing insertion states ============================================================

    // === Recurse, then continue with a deletion afterwards =====================================
    // --- First we put the deletion continuation on the stack...   ------------------------------
    PUSH(DELETE(S.left, S.right, S.pos, S.d, S.edits, /* a=*/1));

    // --- Then we continue with the recursive operation on pos-1 ------------------------------------
    S.pos--;
    S.d--;
    *S.edits++ = 'I';
    goto rec;

    // clang-format off
    // stupid formatter thinks this is the C++ delete...
delete:  // === Processing deletion states ==============================================================
    // clang-format on 
    // If there are no more letters for deletions, or if we are in the first operation, we are done here.
    if (S.a == itr->preproc->alpha.size || S.edits == itr->edits_buf)
    {
        goto pop_next;
    }

    // === Recurse, then continue with a deletion afterwards =====================================
    // --- First we put the deletion continuation on the stack...   ------------------------------
    PUSH(DELETE(S.left, S.right, S.pos, S.d, S.edits, (uint8_t)(S.a + 1)));

    // --- Then we continue with the recursive operation on pos-1 --------------------------------
    S.left = C(S.a) + O(S.a, S.left);
    S.right = C(S.a) + O(S.a, S.right);
    S.d--;
    *S.edits++ = 'D';
    goto rec;

emit:
    if (S.next == S.end)
        goto pop_next; // Done emitting, so continue search

    // Otherwise, emit a hit.
    // Since we return, we have to make sure that we are back in EMIT when we dispatch
    // the next time. For the other transitions, it only maters when we push, but returning
    // works a bit like pushing. If we don't set the state, we will dispatch to a random
    // one when we are called again...
    S.state = EMIT_;
    return (cstr_approx_match){
        .pos = itr->preproc->sa->buf[S.next++], .cigar = itr->cigar_buf};

done: // === If we end up here, there is nothing more to iterate =============================
    return (cstr_approx_match){ .pos = -1, .cigar = ""};
}
#undef S
#undef PUSH
#undef POP

cstr_approx_matcher *cstr_li_durbin_search(cstr_li_durbin_preproc *preproc,
                                           cstr_const_sslice p, long long d)
{
    cstr_approx_matcher *itr = cstr_malloc(sizeof *itr);
    itr->preproc = preproc;
    // We can have p.len match plus d other operations (and '\0')
    size_t edit_len = (size_t)p.len + (size_t)d + 1llu;
    itr->edits_buf = malloc(edit_len); // cstr_alloc_sslice_buf(0, p.len + d);
    // The cigar can at most be twice the length of the edits;
    itr->cigar_buf = malloc(2 * edit_len);

    itr->p_buf = cstr_alloc_sslice(p.len);
    itr->p = CSTR_SLICE_CONST_CAST(*itr->p_buf);
    bool map_ok = cstr_alphabet_map(*itr->p_buf, p, &preproc->alpha);
    itr->stack = new_stack();
    itr->current_state = DONE(); // In case we can't map, we want to be done...

    // If mapping failed, we leave the stack empty. Then the iterator will
    // return no matches, but it will still be in a state where we can free
    // it as per usual.

    if (map_ok)
    {
        // Put a DONE on the stack so we know when to stop, then start with a
        // REC state.
        push_frame(&itr->stack, DONE()); // Marker for when we complete
        itr->current_state = REC(0, preproc->sa->len, p.len - 1, d, itr->edits_buf);
    }

    return itr;
}

void cstr_free_approx_matcher(cstr_approx_matcher *matcher)
{
    free(matcher->p_buf);
    free(matcher->edits_buf);
    free(matcher->cigar_buf);
    free(matcher->stack);
    free(matcher);
}

// Serialisation

static inline size_t c_table_size(size_t sigma)
{
    struct c_table dummy; // To get the size of elements
    return offsetof(struct c_table, cumsum) + sizeof(dummy.cumsum[0]) * sigma;
}

static inline size_t o_table_size(size_t len, size_t sigma)
{
    struct o_table dummy; // To get the size of elements
    return offsetof(struct o_table, table) + sizeof(dummy.table[0]) * sigma * len;
}

void cstr_write_ld_tables(FILE *f, cstr_li_durbin_preproc *tables)
{
    size_t len = (size_t)tables->sa->len;
    fwrite(&len, sizeof len, 1, f);
    fwrite(&tables->alpha, sizeof tables->alpha, 1, f);
    fwrite(tables->sa->buf, len * sizeof tables->sa->buf[0], 1, f);
    fwrite(tables->ctab, c_table_size(tables->alpha.size), 1, f);
    fwrite(tables->otab, o_table_size(len, tables->alpha.size), 1, f);
    fwrite(tables->rotab, o_table_size(len, tables->alpha.size), 1, f);
}

cstr_li_durbin_preproc *cstr_read_ld_tables(FILE *f)
{
    size_t len;
    cstr_li_durbin_preproc *tables = cstr_malloc(sizeof *tables);

    size_t read;
    read = fread(&len, sizeof len, 1, f);
    assert(read == 1);
    read = fread(&tables->alpha, sizeof tables->alpha, 1, f);
    assert(read == 1);
    tables->sa = cstr_alloc_uislice((long long)len);
    read = fread(tables->sa->buf, len * sizeof tables->sa->buf[0], 1, f);
    assert(read == 1);
    tables->ctab = cstr_malloc(c_table_size(tables->alpha.size));
    read = fread(tables->ctab, c_table_size(tables->alpha.size), 1, f);
    assert(read == 1);
    tables->otab = cstr_malloc(o_table_size(len, tables->alpha.size));
    read = fread(tables->otab, o_table_size(len, tables->alpha.size), 1, f);
    assert(read == 1);
    tables->rotab = cstr_malloc(o_table_size(len, tables->alpha.size));
    read = fread(tables->rotab, o_table_size(len, tables->alpha.size), 1, f);
    assert(read == 1);

    return tables;
}

static void print_c_table(struct c_table const *ctab)
{
    printf("[ ");
    for (long long i = 0; i < ctab->sigma; i++)
    {
        printf("%u ", ctab->cumsum[i]);
    }
    printf("]\n");
}

static void print_o_table(struct o_table const *otab)
{
    for (int a = 0; a < otab->sigma; a++)
    {
        for (int i = 0; i < otab->n; i++)
        {
            printf("%lld ", O_RAW(a, i));
        }
        printf("\n");
    }
}

TL_TEST(build_ld_tables)
{
    TL_BEGIN();

    cstr_const_sslice x = CSTR_SLICE_STRING0((const char *)"mississippi");
    cstr_li_durbin_preproc *preproc = cstr_li_durbin_preprocess(x);
    printf("SA: ");
    CSTR_SLICE_PRINT(*preproc->sa);
    printf("\n");
    printf("C: ");
    print_c_table(preproc->ctab);
    printf("\n");
    printf("O:\n");
    print_o_table(preproc->otab);
    printf("\n");
    printf("RO:\n");
    print_o_table(preproc->rotab);

    cstr_free_li_durbin_preproc(preproc);

    TL_END();
}

TL_TEST(ld_iterator)
{
    TL_BEGIN();

    cstr_const_sslice x = CSTR_SLICE_STRING0((const char *)"mississippi");
    cstr_const_sslice p = CSTR_SLICE_STRING((const char *)"ssi");
    cstr_li_durbin_preproc *preproc = cstr_li_durbin_preprocess(x);
    long long d = 1; // FIXME: FOR NOW
    cstr_approx_matcher *itr = cstr_li_durbin_search(preproc, p, d);

    for (cstr_approx_match m = cstr_approx_next_match(itr); m.pos != -1; m = cstr_approx_next_match(itr))
    {
        printf("%lld: %s\n", m.pos, m.cigar);
    }

    cstr_free_approx_matcher(itr);

    TL_END();
}
