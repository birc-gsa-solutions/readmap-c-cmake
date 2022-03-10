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

// MARK: Continuation/closure boiler plate for continuation-passing-style recursion

// We use continuations to avoid exhausting the call-stack. With an optimising
// compiler, all the calls are tail-calls, optimised into jump instructions,
// and the memory usage is all on the explicit stack. The stack contains callback functions
// (continuation) and data for them (closures). When a function cannot do anything more, it should
// call the next continuation (call_continuation). It is a little more primitive than
// full continuation-passing-style, because we only need a stack, so we don't need to
// store continuations in the closures, but can always get the next from the stack
// instead. Thus, when we want to call with a continuation, we instead push the continuation
// to the stack and then call. The continuation will be called when the full processing
// of the called function is done.
union closure;
struct context;
struct stack;
struct stack_frame;

static void push_frame(struct stack **stack, struct stack_frame frame);
static struct stack_frame *pop_frame(struct stack **stack);

// Data used for all search functions.
struct context
{
    cstr_li_durbin_preproc *preproc;
    cstr_sslice *p_buf;      // buffer for alphabet-mapped string
    cstr_const_sslice p;     // slice for p_buf
    cstr_sslice_buf *edits;  // edits so far in the recursion
    char *cigar_buf;         // buffer for the cigar string
    cstr_approx_match match; // used for returning results
    struct stack *stack;     // recursion stack (simple mem management for CPS)
};

// Continuations are really trampoline-thunks that return true if we should continue
// and false if we should not.
typedef bool (*continuation)(struct context *context, union closure *cl);
#define RETURN(POS, CIGAR)                                                      \
    do                                                                          \
    {                                                                           \
        context->match = ((cstr_approx_match){.pos = (POS), .cigar = (CIGAR)}); \
        return false;                                                           \
    } while (0)

#define RETURN_THEN(POS, CIGAR, K)                                 \
    do                                                                          \
    {                                                                           \
        push_frame(&context->stack, K);                                         \
        context->match = ((cstr_approx_match){.pos = (POS), .cigar = (CIGAR)}); \
        return false;                                                           \
    } while (0)

#define CALL(K)            \
    do                                  \
    {                                   \
        push_frame(&context->stack, K); \
        return true;                    \
    } while (0)

#define NEXT() \
    do                           \
    {                            \
        return true;             \
    } while (0)

#define CALL_THEN(CALL, K)    \
    do                                     \
    {                                      \
        push_frame(&context->stack, K);    \
        push_frame(&context->stack, CALL); \
        return true;                       \
    } while (0)

// Construct a continuation from a function and a closure
#define K(CONT, CL) \
    (struct stack_frame) { .k = (CONT), .cl = (CL) }

// clang-format off
struct emit_closure
{
    long long next;
    long long end;
    const char *cigar;
};
#define EMIT_CLOSURE(NEXT, END, CIGAR)   \
    (union closure) { .emit_closure = {  \
        .next = (NEXT), .end = (END),    \
        .cigar = (CIGAR)                 \
    } }

static bool emit(struct context *context, union closure *closure);
#define EMIT(NEXT, END, CIGAR) K(emit, EMIT_CLOSURE(NEXT, END, CIGAR))

struct match_closure
{
    long long left;   // Range where we have matching 
    long long right;  // prefixes.
    long long i;      // Index into p
    long long d;      // Number of edits left
    uint8_t a;        // Current letter we are matching/mismatching
    cstr_sslice_buf_slice edits; // Edit operations so far
};
#define MATCH_CLOSURE(LEFT, RIGHT, I, D, A, EDITS)  \
    (union closure) { .match_closure = {            \
        .left = (LEFT), .right = (RIGHT),           \
        .i = (I), .d = (D), .a = (A),               \
        .edits = (EDITS)                            \
    } }

static bool match(struct context *context, union closure *closure);
#define MATCH(LEFT, RIGHT, I, D, A, EDITS) K(match, MATCH_CLOSURE(LEFT, RIGHT, I, D, A, EDITS))

struct insert_closure
{
    long long left;   // Range where we have matching
    long long right;  // prefixes.
    long long i;      // Index into p
    long long d;      // Number of edits left
    cstr_sslice_buf_slice edits; // Edit operations so far
};
#define INSERT_CLOSURE(LEFT, RIGHT, I, D, EDITS) \
    (union closure) { .insert_closure = {        \
        .left = (LEFT), .right = (RIGHT),        \
        .i = (I), .d = (D),                      \
        .edits = (EDITS)                         \
    } }

static bool insert(struct context *context, union closure *closure);
#define INSERT(LEFT, RIGHT, I, D, EDITS) K(insert, INSERT_CLOSURE(LEFT, RIGHT, I, D, EDITS))

struct delete_closure
{
    long long left;   // Range where we have matching
    long long right;  // prefixes.
    long long i;      // Index into p
    long long d;      // Number of edits left
    uint8_t a;        // Current letter we are "deleting"
    cstr_sslice_buf_slice edits; // Edit operations so far
};
#define DELETE_CLOSURE(LEFT, RIGHT, I, D, A, EDITS)  \
    (union closure) { .delete_closure = {            \
        .left = (LEFT), .right = (RIGHT),            \
        .i = (I), .d = (D), .a = (A),                \
        .edits = (EDITS)                             \
    } }

static bool delete(struct context *context, union closure *closure);
#define DELETE(LEFT, RIGHT, I, D, A, EDITS) K(delete, DELETE_CLOSURE(LEFT, RIGHT, I, D, A, EDITS))

struct rec_search_closure
{
    long long left;   // Range where we have matching 
    long long right;  // prefixes.
    long long i;      // Index into p
    long long d;      // Number of edits left
    cstr_sslice_buf_slice edits; // Edit operations so far
};
#define REC_SEARCH_CLOSURE(LEFT, RIGHT, I, D, EDITS)  \
    (union closure) { .rec_search_closure = {         \
        .left = (LEFT), .right = (RIGHT),             \
        .i = (I), .d = (D),                           \
        .edits = (EDITS)                              \
    } }

static bool rec_search(struct context *context, union closure *closure);
#define REC_SEARCH(LEFT, RIGHT, I, D, EDITS) K(rec_search, REC_SEARCH_CLOSURE(LEFT, RIGHT, I, D, EDITS))


union closure
{
    struct emit_closure emit_closure;
    struct match_closure match_closure;
    struct insert_closure insert_closure;
    struct delete_closure delete_closure;
    struct rec_search_closure rec_search_closure;
};
// clang-format on

struct stack_frame
{
    continuation k;
    union closure cl;
};

// MARK: Stack used for CPS.
// Quite simple implementation, but it suffices for what we need here.

struct stack
{
    size_t size;
    size_t used;
    struct stack_frame frames[];
};

// We only shrink a stack when it is 1/4 used, and then only to 1/2 size,
// so memory we pop off is still available until the next stack action.
#define MIN_STACK_SIZE 128
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

static inline void push_frame(struct stack **stack, struct stack_frame k)
{
    *stack = resize_stack(*stack);
    (*stack)->frames[(*stack)->used++] = k;
}

static inline struct stack_frame *pop_frame(struct stack **stack)
{
    *stack = resize_stack(*stack);
    return &(*stack)->frames[--(*stack)->used];
}

// MARK: actions in the approximative search
static bool emit(struct context *context, union closure *cl)
{
    struct emit_closure *ecl = &cl->emit_closure;
    long long next = ecl->next, end = ecl->end;
    const char *cigar = ecl->cigar;

    if (next < end)
    {
        RETURN_THEN(
            context->preproc->sa->buf[next], cigar,
            EMIT(next + 1, end, cigar));
    }
    else
    {
        NEXT();
    }
}
static bool match(struct context *context, union closure *cl)
{
    struct match_closure *mcl = &cl->match_closure;
    long long left = mcl->left, right = mcl->right, i = mcl->i, d = mcl->d;
    uint8_t a = mcl->a;
    cstr_sslice_buf_slice edits = mcl->edits;

    if (a == context->preproc->alpha.size)
    {
        // No more match operations. Try insertions.
        CALL(INSERT(left, right, i, d, edits));
    }
    else
    {
        struct c_table *ctab = context->preproc->ctab;
        struct o_table *otab = context->preproc->otab;
        long long new_left = C(a) + O(a, left);
        long long new_right = C(a) + O(a, right);
        long long new_d = d - (context->p.buf[i] != a);
        // Recurse, then continue with the next match afterwards
        CALL_THEN(
            REC_SEARCH(new_left, new_right, i - 1, new_d, CSTR_BUF_SLICE_APPEND(edits, 'M')),
            MATCH(left, right, i, d, (uint8_t)(a + 1), edits));
    }
}

static inline bool insert(struct context *context, union closure *cl)
{
    struct insert_closure *icl = &cl->insert_closure;
    long long left = icl->left, right = icl->right, i = icl->i, d = icl->d;
    cstr_sslice_buf_slice edits = icl->edits;

    // Recurse, then continue with deletion afterwards
    CALL_THEN(
        REC_SEARCH(left, right, i - 1, d - 1, CSTR_BUF_SLICE_APPEND(edits, 'I')),
        DELETE(left, right, i, d, (uint8_t)1, edits));
}

static inline bool delete(struct context *context, union closure *cl)
{
    struct delete_closure *dcl = &cl->delete_closure;
    long long left = dcl->left, right = dcl->right, i = dcl->i, d = dcl->d;
    uint8_t a = dcl->a;
    cstr_sslice_buf_slice edits = dcl->edits;

    if (a == context->preproc->alpha.size)
    {
        // No more delete operations. We are done in this path
        // of the exploration
        NEXT();
    }
    if (edits.from == edits.to)
    {
        // If we don't have any edits yet, we don't delete.
        // Those deletions are not interesting.
        NEXT();
    }

    struct c_table *ctab = context->preproc->ctab;
    struct o_table *otab = context->preproc->otab;
    long long new_left = C(a) + O(a, left);
    long long new_right = C(a) + O(a, right);
    // Recurse, then continue afterwards
    CALL_THEN(
        REC_SEARCH(new_left, new_right, i, d - 1, CSTR_BUF_SLICE_APPEND(edits, 'D')),
        DELETE(left, right, i, d, (uint8_t)(a + 1), edits));
}

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

static bool rec_search(struct context *context, union closure *cl)
{
    struct rec_search_closure *rscl = &cl->rec_search_closure;
    long long left = rscl->left, right = rscl->right, i = rscl->i, d = rscl->d;
    cstr_sslice_buf_slice edits = rscl->edits;
    
    if (left >= right || d < 0)
    {
        // Nothing to be found here, try the next continuation
        NEXT();
    }
    if (i < 0)
    {
        // We have a match, so emit it
        context->cigar_buf = cstr_realloc(context->cigar_buf, (size_t)(2 * context->edits->cap + 1));
        edits_to_cigar(context->cigar_buf, CSTR_SLICE_CONST_CAST(CSTR_BUF_SLICE_DEREF(edits)));
        CALL(EMIT(left, right, context->cigar_buf));
    }

    // Otherwise, continue the recursion with the next match operation.
    // The match operation starts with a == 1 because we don't want
    // the sentinel.
    CALL(MATCH(left, right, i, d, /* a = */ 1, edits));
}

static inline bool bounce(struct context *context)
{
    if (context->stack->used == 0)
    {
        // NO MORE MATCHES
        RETURN(-1, "");
    }
    struct stack_frame *sf = pop_frame(&context->stack);
    return sf->k(context, &sf->cl);
}

static cstr_approx_match trampoline(struct context *context)
{
    // Bounce until there is no more bouncing
    while (bounce(context))
        ;
    // and then return the result from context
    return context->match;
}

// If I implement other approximative matcheres, then
// this should be the Li & Durbin specific, and we need a
// vtab as for the exact matcher to distinguish between
// them.
struct cstr_approx_matcher
{
    cstr_sslice *p_buf;
    struct context context;
};

cstr_approx_matcher *cstr_li_durbin_search(cstr_li_durbin_preproc *preproc,
                                           cstr_const_sslice p, long long d)
{
    cstr_approx_matcher *itr = cstr_malloc(sizeof *itr);
    itr->context.preproc = preproc;
    itr->context.edits = cstr_alloc_sslice_buf(0, p.len + d);
    itr->context.cigar_buf = 0; // We (re)alloc when we need it

    itr->p_buf = cstr_alloc_sslice(p.len);
    itr->context.p = CSTR_SLICE_CONST_CAST(*itr->p_buf);
    bool map_ok = cstr_alphabet_map(*itr->p_buf, p, &preproc->alpha);
    itr->context.stack = new_stack();

    if (map_ok)
    {
        // If mapping failed, we leave the stack empty. Then the iterator will
        // return no matches, but it will still be in a state where we can free
        // it as per usual.
        cstr_sslice_buf_slice edits = {.buf = &itr->context.edits, .from = 0, .to = 0};
        push_frame(&itr->context.stack,
                   REC_SEARCH(0, preproc->sa->len, p.len - 1, d, edits));
    }

    return itr;
}

void cstr_free_approx_matcher(cstr_approx_matcher *matcher)
{
    free(matcher->p_buf);
    free(matcher->context.edits);
    free(matcher->context.cigar_buf);
    free(matcher->context.stack);
    free(matcher);
}

cstr_approx_match cstr_approx_next_match(cstr_approx_matcher *matcher)
{
    return trampoline(&matcher->context);
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

#ifdef GEN_UNIT_TESTS // unit testing of static functions...

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

#endif // GEN_UNIT_TESTS
