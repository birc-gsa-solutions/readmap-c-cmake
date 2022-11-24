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

static inline long long scan_next(long long i, cstr_const_sslice edits)
{
    long long j = i;
    for (; j >= 0; j--)
    {
        if (edits.buf[i] != edits.buf[j])
        {
            return j;
        }
    }
    // If we get here, we scanned the last segment and they are all equal,
    // so we should return -1 as past the last (in reverse order)
    return -1;
}

static void edits_to_cigar(char *cigar_buf, cstr_const_sslice edits)
{
    size_t space_left = 2 * (size_t)edits.len + 1; // We allocated at least this much
    for (long long i = edits.len - 1, j = scan_next(i, edits); i >= 0; i = j, j = scan_next(i, edits))
    {
        int used = snprintf(cigar_buf, space_left, "%d%c", (int)(i - j), (char)edits.buf[i]);
        cigar_buf = cigar_buf + used;
        space_left -= (size_t)used;
    }
    *cigar_buf = '\0';
}

// MARK: State machine for searching
struct done_state
{
    // When we are done, we don't do anything, so we don't
    // need to remember any data...
};
struct emit_state
{
    long long next;    // Next index in sa to emit
    long long end;     // Where the hits end
    const char *cigar; // The cigar that matches the match
};
struct rec_state
{
    long long left;              // Range where we have matching
    long long right;             // prefixes.
    long long pos;               // Index into p
    long long d;                 // Number of edits left
    cstr_sslice_buf_slice edits; // Edit operations so far
};
struct match_state
{
    long long left;              // Range where we have matching
    long long right;             // prefixes.
    long long pos;               // Index into p
    long long d;                 // Number of edits left
    cstr_sslice_buf_slice edits; // Edit operations so far
    uint8_t a;                   // Current letter we are matching/mismatching
};
struct insert_state
{
    long long left;              // Range where we have matching
    long long right;             // prefixes.
    long long pos;               // Index into p
    long long d;                 // Number of edits left
    cstr_sslice_buf_slice edits; // Edit operations so far
};
struct delete_state
{
    long long left;              // Range where we have matching
    long long right;             // prefixes.
    long long pos;               // Index into p
    long long d;                 // Number of edits left
    cstr_sslice_buf_slice edits; // Edit operations so far
    uint8_t a;                   // Current letter we are "deleting"
};

struct state;
typedef struct state state;
typedef bool (*transition_fn)(cstr_approx_matcher *itr, state *s, cstr_approx_match *m);
struct state
{
    union
    {
        struct done_state done;     // Done is just a marker to signal termination
        struct emit_state emit;     // Emit outputs matches one at a time
        struct rec_state rec;       // Rec starts the recursion (M->I->D) at a given state
        struct match_state match;   // Match handles a single match
        struct insert_state insert; // Insert tries an insertion
        struct delete_state delete; // and Delete handles a deletion
    };
    transition_fn fn;
};

static inline bool next(cstr_approx_matcher *itr, state *s, cstr_approx_match *m)
{
    // Call the state's transition function.
    // We don't need a swithc when we have function calls...
    // (Assuming that tail call optimisation works, of course).
    // The function basically works like a trampoline in
    // continuation-passing-style programming except that
    // we use an explicit stack.
    return s->fn(itr, s, m);
}

// clang-format off
static bool done_transition  (cstr_approx_matcher *itr, state *s, cstr_approx_match *m);
static bool emit_transition  (cstr_approx_matcher *itr, state *s, cstr_approx_match *m);
static bool rec_transition   (cstr_approx_matcher *itr, state *s, cstr_approx_match *m);
static bool match_transition (cstr_approx_matcher *itr, state *s, cstr_approx_match *m);
static bool insert_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m);
static bool delete_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m);
// clang-format on

#define DONE() ((state){.fn = done_transition})
#define EMIT(NEXT, END, CIGAR)                                   \
    ((state){.emit = {.next = NEXT, .end = END, .cigar = CIGAR}, \
             .fn = (transition_fn)emit_transition})
#define REC(LEFT, RIGHT, POS, D, EDITS)                                                 \
    ((state){.rec = {.left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS}, \
             .fn = (transition_fn)rec_transition})
#define MATCH(LEFT, RIGHT, POS, D, EDITS, A)                                                      \
    ((state){.match = {.left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS, .a = A}, \
             .fn = (transition_fn)match_transition})
#define INSERT(LEFT, RIGHT, POS, D, EDITS)                                                 \
    ((state){.insert = {.left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS}, \
             .fn = (transition_fn)insert_transition})
#define DELETE(LEFT, RIGHT, POS, D, EDITS, A)                                                     \
    ((state){.match = {.left = LEFT, .right = RIGHT, .pos = POS, .d = D, .edits = EDITS, .a = A}, \
             .fn = (transition_fn)delete_transition})

#define PUSH(STATE) push_frame(&itr->stack, STATE)
#define POP() *(pop_frame(&itr->stack))
#define RUN_STACK_TOP() \
    do                  \
    {                   \
        *s = POP();     \
        return false;   \
    } while (0)
#define RUN_NEXT_STATE(STATE) \
    do                        \
    {                         \
        *s = STATE;           \
        return false;         \
    } while (0)

// MARK: Stack used for CPS.
// Quite simple implementation, but it suffices for what we need here.

// FIXME: I should be able to work out the max stack depth a priori so I don't have to
// waste time on resizing...
struct stack
{
    size_t size;
    size_t used;
    state frames[];
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
    if (stack->size / 4 < stack->used && stack->used < stack->size)
    {
        return stack; // No resize necessary
    }
    if (stack->used < MIN_STACK_SIZE)
    {
        // Don't shrink below this size
        assert(stack->size > stack->used);
        return stack;
    }
    if (stack->used == stack->size)
    {
        stack->size *= 2;
    }
    else
    {
        assert(stack->used <= stack->size / 4);
        stack->size /= 2;
        stack->size = (stack->size < MIN_STACK_SIZE) ? MIN_STACK_SIZE : stack->size;
    }
    return cstr_realloc_header_array(stack,
                                     offsetof(struct stack, frames),
                                     sizeof stack->frames[0],
                                     stack->size);
}

static inline void push_frame(struct stack **stack, state state)
{
    *stack = resize_stack(*stack);
    (*stack)->frames[(*stack)->used++] = state;
}

static inline state *pop_frame(struct stack **stack)
{
    *stack = resize_stack(*stack);
    return &(*stack)->frames[--(*stack)->used];
}

// MARK: iterator for searches
// If I implement other approximative matcheres, then
// this should be the Li & Durbin specific, and we need a
// vtab as for the exact matcher to distinguish between
// them.
struct cstr_approx_matcher
{
    state current_state;
    cstr_sslice *p_buf;
    cstr_const_sslice p;
    cstr_li_durbin_preproc *preproc;
    cstr_sslice_buf *edits_buf;
    char *cigar_buf;
    struct stack *stack;
};

// MARK: State transitions
static bool done_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m)
{
    // When we are done, tell the matcher and return true so iterator returns
    *m = (cstr_approx_match){.pos = -1, .cigar = ""};
    return true;
}

static bool emit_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m)
{
    if (s->emit.next != s->emit.end)
    {
        // Emit a match if there are more...
        *m = (cstr_approx_match){.pos = s->emit.next++, .cigar = s->emit.cigar};
        return true;
    }
    // Otherwise, the next state is at the top of the stack
    *s = POP();
    return false; // let the next state take over
}

static bool rec_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m)
{
    if (s->rec.left == s->rec.right || s->rec.d < 0)
    {
        // No matches here, so continue with the next continuation...
        RUN_STACK_TOP();
    }

    if (s->rec.pos < 0)
    {
        // We have a match!
        itr->cigar_buf = cstr_realloc(itr->cigar_buf, (size_t)(2 * itr->edits_buf->cap + 1));
        edits_to_cigar(itr->cigar_buf, CSTR_SLICE_CONST_CAST(CSTR_BUF_SLICE_DEREF(s->rec.edits)));
        RUN_NEXT_STATE(EMIT(s->rec.left, s->rec.right, itr->cigar_buf));
    }

    // Otherwise, continue the recursion with the next match operation.
    // The match operation starts with a == 1 because we don't want
    // the sentinel.
    RUN_NEXT_STATE(MATCH(s->rec.left, s->rec.right, s->rec.pos, s->rec.d, s->rec.edits, /* a = */ 1));
}

static bool match_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m)
{
    // We pull these into the name space for the O() and C() macros.
    struct c_table *ctab = itr->preproc->ctab;
    struct o_table *otab = itr->preproc->otab;

    if (s->match.a == itr->preproc->alpha.size)
    {
        // We are through the alphabet, so there is nothing more to match.
        // Try inserting instead.
        RUN_NEXT_STATE(INSERT(s->match.left, s->match.right, s->match.pos, s->match.d, s->match.edits));
    }

    // === Recurse, then continue with the next match afterwards =====================================
    // --- First we put the match continuation on the stack...   -------------------------------------
    PUSH(MATCH(s->match.left, s->match.right, s->match.pos, s->match.d, s->match.edits, s->match.a + 1));
    // --- Then we continue with the recursive operation on pos-1 ------------------------------------
    long long new_left = C(s->match.a) + O(s->match.a, s->match.left);
    long long new_right = C(s->match.a) + O(s->match.a, s->match.right);
    long long new_d = s->match.d - (itr->p.buf[s->match.pos] != s->match.a);
    RUN_NEXT_STATE(REC(new_left, new_right, s->match.pos - 1, new_d,
                       CSTR_BUF_SLICE_APPEND(s->match.edits, 'M')));
}

static bool insert_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m)
{
    // === Recurse, then continue with a deletion afterwards =====================================
    // --- First we put the deletion continuation on the stack...   ------------------------------
    PUSH(DELETE(s->insert.left, s->insert.right, s->insert.pos, s->insert.d, s->insert.edits, /* a=*/1));
    // --- Then we continue with the recursive operation on pos-1 ------------------------------------
    RUN_NEXT_STATE(REC(s->insert.left, s->insert.right, s->insert.pos - 1, s->insert.d - 1,
                       CSTR_BUF_SLICE_APPEND(s->insert.edits, 'I')));
}

static bool delete_transition(cstr_approx_matcher *itr, state *s, cstr_approx_match *m)
{
    unsigned char a = s->delete.a;
    if (a == itr->preproc->alpha.size)
    {
        // No more delete operations. We are done in this path
        // of the exploration
        RUN_STACK_TOP();
    }
    if (s->delete.edits.from == s->delete.edits.to)
    {
        // If we don't have any edits yet, we don't delete.
        // Those deletions are not interesting.
        RUN_STACK_TOP();
    }

    struct c_table *ctab = itr->preproc->ctab;
    struct o_table *otab = itr->preproc->otab;
    long long left = s->delete.left, right = s->delete.right;
    long long new_left = C(a) + O(a, left);
    long long new_right = C(a) + O(a, right);
    // === Recurse, then continue with a deletion afterwards =====================================
    // --- First we put the deletion continuation on the stack...   ------------------------------
    PUSH(DELETE(left, right, s->insert.pos, s->insert.d - 1, s->insert.edits, a + 1));
    // --- Then we continue with the recursive operation on pos-1 ------------------------------------
    RUN_NEXT_STATE(REC(new_left, new_right, s->insert.pos, s->insert.d - 1,
                       CSTR_BUF_SLICE_APPEND(s->insert.edits, 'D')));
}

cstr_approx_matcher *cstr_li_durbin_search(cstr_li_durbin_preproc *preproc,
                                           cstr_const_sslice p, long long d)
{
    cstr_approx_matcher *itr = cstr_malloc(sizeof *itr);
    itr->preproc = preproc;
    itr->edits_buf = cstr_alloc_sslice_buf(0, p.len + d);
    itr->cigar_buf = 0; // We (re)alloc when we need it

    itr->p_buf = cstr_alloc_sslice(p.len);
    itr->p = CSTR_SLICE_CONST_CAST(*itr->p_buf);
    bool map_ok = cstr_alphabet_map(*itr->p_buf, p, &preproc->alpha);
    itr->stack = new_stack();
    itr->current_state = DONE(); // In case we can't map, we want to be done...

    if (map_ok)
    {
        // If mapping failed, we leave the stack empty. Then the iterator will
        // return no matches, but it will still be in a state where we can free
        // it as per usual.
        cstr_sslice_buf_slice edits = {.buf = &itr->edits_buf, .from = 0, .to = 0};
        push_frame(&itr->stack, DONE()); // Marker for when we complete
        itr->current_state = REC(0, preproc->sa->len, p.len - 1, d, edits);
    }

    return itr;
}

cstr_approx_match cstr_approx_next_match(cstr_approx_matcher *matcher)
{
    cstr_approx_match match = {.pos = -1, .cigar = ""};
    while (next(matcher, &matcher->current_state, &match)) // FIXME: DO I NEED TO LOOP HERE OR CAN I CALL DIRECTLY IN TRANSITIONS?
        /* Run like a trampoline until we have to return something. */;

#if 0
    // We pull these into the name space for the O() and C() macros.
    struct c_table *ctab = matcher->preproc->ctab;
    struct o_table *otab = matcher->preproc->otab;

    cstr_approx_match match = {.pos = -1, .cigar = ""};
    while (matcher->stack->used > 0)
    {
        struct state *state = pop_frame(&matcher->stack);
        switch (state->op)
        {


        case M:
        {

        break;

        case I:
        {

        }
        break;

        case D:
        {
            if (state->a == matcher->preproc->alpha.size)
            {
                // No more delete operations. We are done in this path
                // of the exploration
                continue;
            }

            if (state->edits.from == state->edits.to)
            {
                // If we don't have any edits yet, we don't delete.
                // Those deletions are not interesting.
                continue;
            }
        }
        break;
        }
    }
#endif

    return match;
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
