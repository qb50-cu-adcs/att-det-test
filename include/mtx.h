#define MAX_ROWS 8
#define MAX_COLS 4

struct mtx_matrix{
    int rows;
    int cols;
    float data[MAX_ROWS*MAX_COLS];
};

void mtx_create(int rows, int cols, float *data, struct mtx_matrix* mtx_out);
void mtx_create_ones(int rows, int cols, struct mtx_matrix* mtx_out);
float mtx_get(int row, int col, struct mtx_matrix* mtx_a);
void mtx_set(int row, int col, struct mtx_matrix* mtx_a, float data);
void mtx_print(struct mtx_matrix* mtx_a);
int mtx_sum(struct mtx_matrix* mtx_a, struct mtx_matrix* mtx_b, struct mtx_matrix* mtx_out);
int mtx_scale(struct mtx_matrix* mtx_a, float scaler, struct mtx_matrix* mtx_out);
int mtx_mult(struct mtx_matrix* mtx_a, struct mtx_matrix* mtx_b, struct mtx_matrix* mtx_out);
int mtx_trans(struct mtx_matrix* mtx_a, struct mtx_matrix* mtx_out);
float mtx_det(struct mtx_matrix* mtx_a);
float mtx_norm(struct mtx_matrix* mtx_a);
