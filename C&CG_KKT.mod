##Set
set SUPPLY;
set DEMAND;

##Param
param f   {SUPPLY};
param a   {SUPPLY};
param k   {SUPPLY};
param c   {SUPPLY,DEMAND};
param d_L {DEMAND}; #basic demand
param d_T {DEMAND}; #maximal deviation

##Other param
param N;  #number of iterations
param LB; #lower bounds
param UB; #upper bounds
param M;  #big-M

##Master Problem Param and Var
param d_star  {DEMAND,1..N};  # from Sub problem
var y        {SUPPLY} binary;
var z        {SUPPLY}>=0;
var ¦Ç;
#column generation
var x_gen   {SUPPLY,DEMAND,1..N}>=0;  

##Sub Problem Param and Var
param z_star   {SUPPLY};  #from Master problem
var x        {SUPPLY,DEMAND};
var g        {DEMAND}>=0;
var d        {DEMAND};
#dual var
var ¦È        {DEMAND}>=0;
var ¦Ð        {SUPPLY}>=0;
#linearization aux var
var v        {SUPPLY}binary;
var w        {DEMAND}binary;
var h        {SUPPLY,DEMAND}binary;

##Master problem
minimize MP:  
            sum{i in SUPPLY} (f[i]*y[i] + a[i]*z[i]) + ¦Ç;

subject to C1  {i in SUPPLY}:   
                               z[i] <= k[i]*y[i];
subject to C2:   
                               sum{i in SUPPLY} z[i] >= 772;  #772 == sum of basic demand + 1.8*maximal deviation
#constraint generation
subject to gen_C1  {n in 1..N}: 
                               ¦Ç >= sum{i in SUPPLY,j in DEMAND} c[i,j]*x_gen[i,j,n];
subject to gen_C2  {i in SUPPLY,n in 1..N}:
                               sum{j in DEMAND} x_gen[i,j,n] <= z[i];
subject to gen_C3  {j in DEMAND,n in 1..N}: 
                               sum{i in SUPPLY} x_gen[i,j,n] >= d_star[j,n];

##Sub problem
maximize SP:  
            sum{i in SUPPLY,j in DEMAND} c[i,j]*x[i,j];

subject to C3  {i in SUPPLY}: 
                             sum{j in DEMAND} x[i,j] <= z_star[i];
subject to C4  {j in DEMAND}: 
                             sum{i in SUPPLY} x[i,j] >= d[j];
subject to C5  {i in SUPPLY,j in DEMAND}:
                             ¦È[j] - ¦Ð[i] <= c[i,j];
subject to C6   {j in DEMAND}: 
                             d[j]=d_L[j] + d_T[j]*g[j];
subject to C7:  
                             sum{j in DEMAND:j!='D3'} g[j] <= 1.2;
subject to C8:  
                             sum{j in DEMAND} g[j] <= 1.8;
subject to C9   {j in DEMAND}: 
                             g[j] <= 1;
subject to C10  {i in SUPPLY}: 
                             ¦Ð[i] <= M*v[i];
subject to C11  {i in SUPPLY}: 
                             z_star[i] - sum{j in DEMAND} x[i,j] <= M*(1-v[i]);
subject to C12  {j in DEMAND}:
                             ¦È[j] <= M*w[j];
subject to C13  {j in DEMAND}: 
                             sum{i in SUPPLY} x[i,j] - d[j] <= M*(1-w[j]);
subject to C14  {i in SUPPLY,j in DEMAND}: 
                             x[i,j] <= M*h[i,j];
subject to C15  {i in SUPPLY,j in DEMAND}: 
                             c[i,j] - ¦È[j] + ¦Ð[i] <= M*(1-h[i,j]);



