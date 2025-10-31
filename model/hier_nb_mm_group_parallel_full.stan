functions {
  // 并行似然切片函数（reduce_sum 标准签名：y_slice, start, end, shared...）
  real nb_llk_slice(array[] int y_slice, int start, int end,
                    array[] int store_idx,
                    array[] int mall_of_store,
                    array[] int row_ptr, array[] int col_idx, vector w,
                    vector mall, vector a_store, vector u,
                    real alpha0, vector beta, matrix X, real phi) {
    int P_local = cols(X);
    real lp = 0;
    for (i in 1:size(y_slice)) {
      int n = start + i - 1;
      int s = store_idx[n];
      real eta = alpha0 + mall[ mall_of_store[s] ] + a_store[s];
      if (P_local > 0) eta += dot_product(row(X, n), beta);
      real z = 0;
      for (k in row_ptr[n]:(row_ptr[n+1]-1)) {
        z += w[k] * u[ col_idx[k] ];
      }
      lp += neg_binomial_2_log_lpmf(y_slice[i] | eta + z, phi);
    }
    return lp;
  }
}

data {
  int<lower=1> N;                       // 店-日观测数
  array[N] int<lower=0> y;              // 响应：visit

  int<lower=1> S;                       // 店铺数（全局唯一：mall::store）
  array[N] int<lower=1, upper=S> store_idx;   // 行对应的店索引

  int<lower=1> M;                       // 商场数
  array[S] int<lower=1, upper=M> mall_of_store; // 每个店属于哪个商场

  // 影响者（列）— CSR 稀疏设计（infl 全局唯一：mall::infl）
  int<lower=1> J;                       // 影响者总数
  int<lower=1> K;                       // 非零元素数
  array[N+1] int<lower=1> row_ptr;     // 行指针（1-based, 长度 N+1）
  array[K] int<lower=1, upper=J> col_idx;     // 列索引（1..J）
  vector[K] w;                          // 非零权重（如 log1p(shared_visits)）

  array[J] int<lower=0, upper=1> group_type;   // 影响者分组类型（0/1）

  // 控制变量（可为空矩阵）
  int<lower=0> P;                       // 控制变量列数
  matrix[N, P] X;                       // 设计矩阵（建议已标准化/去均值）
}

parameters {
  real alpha0;                          // 全局截距

  // mall / store 两层随机截距（非中心化）
  vector[M] mall_raw;                   // ~ N(0,1)
  real<lower=0> sigma_mall;

  vector[S] store_raw;                  // ~ N(0,1)
  real<lower=0> sigma_store;

  // 影响者层随机效应（按组不同均值和方差）
  vector[J] u_raw;                      // ~ N(0,1)
  real mu_group0;                       // group_type=0 组均值
  real mu_group1;                       // group_type=1 组均值
  real<lower=0> sigma_group0;           // group_type=0 组 SD
  real<lower=0> sigma_group1;           // group_type=1 组 SD

  // 控制变量系数
  vector[P] beta;

  // 负二项过度分散
  real<lower=0> phi;
}

transformed parameters {
  vector[M] mall   = sigma_mall  * mall_raw;
  vector[S] a_store = sigma_store * store_raw;

  // 组别均值 + 方差（组内收缩）
  vector[J] u;
  for (j in 1:J) {
    if (group_type[j] == 1) {
      u[j] = mu_group1 + sigma_group1 * u_raw[j];
    } else {
      u[j] = mu_group0 + sigma_group0 * u_raw[j];
    }
  }
}

model {
  // 先验（与标准化后的 X / w 相适配，更稳健）
  alpha0       ~ normal(0, 5);

  mall_raw     ~ normal(0, 1);
  sigma_mall   ~ normal(0, 1);

  store_raw    ~ normal(0, 1);
  sigma_store  ~ normal(0, 1);

  u_raw        ~ normal(0, 1);
  mu_group0    ~ normal(0, 1);
  mu_group1    ~ normal(0, 1);
  sigma_group0 ~ student_t(3, 0, 1);  // 厚尾先验
  sigma_group1 ~ student_t(3, 0, 1);  // 厚尾先验

  if (P > 0) beta ~ normal(0, 1);

  phi ~ exponential(1);

  // 并行似然累加（reduce_sum，grainsize=1000）
  target += reduce_sum(nb_llk_slice, y, 1000,
                       store_idx, mall_of_store,
                       row_ptr, col_idx, w,
                       mall, a_store, u,
                       alpha0, beta, X, phi);
}

generated quantities {
  // 便于诊断的线性预测（log 均值）
  vector[N] log_mu_hat;
  for (n in 1:N) {
    int s = store_idx[n];
    real eta = alpha0 + mall[ mall_of_store[s] ] + a_store[s];
    if (cols(X) > 0) eta += dot_product(row(X, n), beta);
    real z = 0;
    for (k in row_ptr[n]:(row_ptr[n+1]-1)) {
      z += w[k] * u[ col_idx[k] ];
    }
    log_mu_hat[n] = eta + z;
  }
}
