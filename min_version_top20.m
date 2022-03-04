clc;
clear;

%% extract data

%load ('GenoData.mat'); % 0-1 genotype information for 200  individuals
%load ('recomData.mat'); % recombination frequency matrix used in "cross" function for producing progeny

% load('GenotypeMulti.mat');
% SelectedGenoData = GenotypeMatrix(1:10000,:,1:200);
% save('GenoData.mat', 'SelectedGenoData');
load('GenoData.mat');
GenoData = SelectedGenoData(1:1000, :,:);
load('recombinationfreq.mat');
recomData = RecombRate(1:999, :);

% beta_height = randfixedsum(1000,1,5000,1,100); % additive effect vector for "height"
% beta_diameter = randfixedsum(1000,1,5000, 1, 100); % additive effect vector for "diameter"
% 
% save('beta_height.mat','beta_height');
% save('beta_diameter.mat','beta_diameter');

load('effects.mat');
beta_height = Height(1:1000,:);
beta_diameter = diameter(1:1000,:);

% in the following loop, we create some weight sets, in each weight set,
% first element is defined for height, second for diameter and third for
% disease resistancy. All the weights are between 0 and 1
r={};
a = 0:0.01:0.999;
a(1) = 0.0001;

b = 1-a;

for k=1:length(a)
    weight = [a(k) b(k)];
    r{k} = weight;
end
% for k =1:100
% r{k} = randfixedsum(3, 1, 1, 0.01, 1)';
% end
% 
% save('/users/mahsa/Box/MT-GS/Codes/DistanceOptimal/r_3D_random.mat', 'r');
fig_id = 1;


linear_pr_pheno = {}; % stores the phenotype information in each generation for linear function
nonlinear_pr_pheno = {}; % stores the phenotype information in each generation for linear function
minh_linear = {};
maxh_linear = {};
minh_nonlinear = {};
maxh_nonlinear = {};

tic

genome_gen_linear ={}; % store geneotype information in each generation for linear function
genome_gen_nonlinear = {}; % store geneotype information in each generation for nonlinear function


pheno_all_linear = {}; % stores the phenotype information with the index of each individual in each generation for linear function
pheno_all_nonlinear = {}; % stores the phenotype information with the index of each individual in each generation for linear function


%% initializing the first generation


genome_gen_linear{1} = GenoData(:,:,:);
genome_gen_nonlinear{1} = GenoData(:,:,:);
parents_linear = {};
parents_nonlinear = {};

% select 20 random parent for each function to start the methahuristic
%algorithm - both functions start with same parents
rnd_parents = randperm(200);
parents_nonlinear{1} = rnd_parents(1:20);
parents_linear{1} = parents_nonlinear{1};



    
    
    point_nonlinear = {};
    
    point_linear = [];
    point_nonlinear = [];
    
    all_points_ri_linear = [];
    all_points_ri_nonlinear = [];
    
    for ri = 1:length(r)
        rep_ri_t1 = [];
        rep_ri_t2 = [];
        rep_ri_t3 = [];
        rep_ri_t1_linear =[];
        rep_ri_t2_linear =[];
        rep_ri_t3_linear =[];
        
        
        for rep = 1:10
            
            
            
            for gen=1:5
                
                %% calculating 3 traits for both linear and nonlinear functions
                [height_pheno_linear, height_pheno_nonlinear] = height_calc(beta_height,genome_gen_linear,genome_gen_nonlinear, gen );
                [diameter_pheno_linear,diameter_pheno_nonlinear ] = diameter_calc(beta_diameter,genome_gen_linear,genome_gen_nonlinear,GenoData, gen );
                 [disease_pheno_linear, disease_pheno_nonlinear] = disease_calc(genome_gen_linear,genome_gen_nonlinear, gen );
                
                
                pheno_all_linear{gen} = [ 1:1:200;height_pheno_linear; diameter_pheno_linear]; % first row defines the index of parents
                pheno_all_nonlinear{gen} = [1:1:200;height_pheno_nonlinear; diameter_pheno_nonlinear];
                
                
                
                % calculating the maximum and minimum of height for normalizing it
%                 min(height_pheno_nonlinear)
%                 max(height_pheno_nonlinear)
                minh_linear{1} = -1;
                maxh_linear{1} = 1;
       
                minh_nonlinear{1} = minh_linear{1};
                maxh_nonlinear{1} = maxh_linear{1};

                mind_linear{1} = -1;
                maxd_linear{1} = 1;
%                 
                mind_nonlinear{1} = mind_linear{1};
                maxd_nonlinear{1} = maxd_linear{1};
                
                
                %%   find 20 best parents for each function out of the metahuristic algorithm
                
                final_parents_nonlinear = find_parents_nonlinear(pheno_all_nonlinear{gen}, gen,minh_nonlinear{1}, maxh_nonlinear{1}, r, ri);
                final_parents_linear = find_parents_linear(pheno_all_linear{gen}, gen, minh_linear{1}, maxh_linear{1},r, ri);

                

                % get the phenotype of current generation
                linear_pr_pheno{gen} = pheno_all_linear{gen}(2:end,:);
                nonlinear_pr_pheno{gen} = pheno_all_nonlinear{gen}(2:end,:);
                
                
                
%                 [minor_allele_linear, minor_allele_nonlinear] = minor_allele_calc(genome_gen_linear, genome_gen_nonlinear, gen);
                
                
            

                
                % creating the next generation
                parents_linear{gen+1} = final_parents_linear;
                parents_nonlinear{gen+1} = final_parents_nonlinear;
                
                %create 5 crosses out of those 20 selected parents and produce 40
                %progeny from each cross
                [pop_linear, pop_nonlinear] = cross(genome_gen_linear,genome_gen_nonlinear, final_parents_linear, final_parents_nonlinear, recomData, gen);
                
                % store the genotype information of the next generation
                genome_gen_linear{gen+1} = pop_linear;
                genome_gen_nonlinear{gen+1} = pop_nonlinear;
                
            end
            

            
            rep_ri_t1 = [rep_ri_t1 , ((nonlinear_pr_pheno{gen}(1,:)-minh_nonlinear{1})/(maxh_nonlinear{1}-minh_nonlinear{1}))];
            rep_ri_t2 = [rep_ri_t2 ,   (nonlinear_pr_pheno{gen}(2,:))];
%             rep_ri_t3 = [rep_ri_t3 ,  mean(nonlinear_pr_pheno{gen}(2,:))];


            rep_ri_t1_linear = [rep_ri_t1_linear , ((linear_pr_pheno{gen}(1,:)-minh_linear{1})/(maxh_linear{1}-minh_linear{1}))];
            rep_ri_t2_linear = [rep_ri_t2_linear , (linear_pr_pheno{gen}(2,:))];
%             rep_ri_t3_linear = [rep_ri_t3_linear ,  mean(linear_pr_pheno{gen}(2,:))];
  
        end
        all_points_ri_linear = [all_points_ri_linear; [(rep_ri_t1_linear'), (rep_ri_t2_linear')] ];
        all_points_ri_nonlinear = [all_points_ri_nonlinear; [(rep_ri_t1'), (rep_ri_t2')]];
        
        
%         point_nonlinear{ri} = [(rep_ri_t1'), (rep_ri_t2')];
%         point_linear{ri} = [(rep_ri_t1_linear'), (rep_ri_t2_linear')];
        ri
    end
%     fname = '/users/mahsa/Box/MT-LATB/Figures/repetitions';
    format short
    save('/users/mahsa/Box/MT-GS/Codes/DistanceOptimal/matrix_rep_min_height_diameter.mat', 'all_points_ri_nonlinear');
    save('/users/mahsa/Box/MT-GS/Codes/DistanceOptimal/matrix_rep_linear_height_diameter.mat', 'all_points_ri_linear');

%     save('C:\Users\famini\Box\MT-GS\Codes\DistanceOptimal\matrix_rep_min_height_diameter.mat', 'all_points_ri_nonlinear');
%     save('C:\Users\famini\Box\MT-GS\Codes\DistanceOptimal\matrix_rep_linear_height_diameter.mat', 'all_points_ri_linear');




%% functions
% height calculation function
function [height_pheno_linear, height_pheno_nonlinear] = height_calc(beta_height,genome_gen_linear,genome_gen_nonlinear, gen )
[m,~,p] = size(genome_gen_linear{gen});

% height is calculated by \beta * G
height_pheno_linear = beta_height' * reshape(sum(genome_gen_linear{gen},2),[m,p]);


height_pheno_nonlinear = beta_height' * reshape(sum(genome_gen_nonlinear{gen},2),[m,p]);

end

% diameter calculation function
function [diameter_pheno_linear,diameter_pheno_nonlinear ] = diameter_calc(beta_diameter,genome_gen_linear,genome_gen_nonlinear, GenoData, gen )

[m,~,p] = size(genome_gen_linear{gen});

% at first diameter is calculated based on \beta *G
diameter_pheno_linear = beta_diameter' *reshape(sum(genome_gen_linear{gen},2),[m,p]);
diameter_pheno_nonlinear = beta_diameter' *reshape(sum(genome_gen_nonlinear{gen},2),[m,p]);

% calculate the 15 %percentile of the diameter of initial population as the lower bound and 75% percentile as
% the upper bound
lowerBound =prctile( beta_diameter' *reshape(sum(GenoData,2),[m,200]),45);
upperBound =prctile( beta_diameter' *reshape(sum(GenoData,2),[m,200]), 75);
% 
% 
% % if lower bound < diameter of each individual < upper bound --->
% % diameter = 1, else diameter = 0
for x=1:length(diameter_pheno_linear)
    if (diameter_pheno_linear(x) >= lowerBound && diameter_pheno_linear(x) <= upperBound)
        diameter_pheno_linear(x) = 1;
    else
        diameter_pheno_linear(x) = 0;
    end
end


for x=1:length(diameter_pheno_nonlinear)
    if (diameter_pheno_nonlinear(x) >= lowerBound && diameter_pheno_nonlinear(x) <= upperBound)
        diameter_pheno_nonlinear(x) = 1;
    else
        diameter_pheno_nonlinear(x) = 0;
    end
end

end


% disease calculation function
function [disease_pheno_linear, disease_pheno_nonlinear] = disease_calc(genome_gen_linear,genome_gen_nonlinear,gen )

% if the sum of the left and right gametes of an individual in location 8 >1 and in
%location 780 > 1 ----> disease = 1, else disease = 0

disease_pheno_linear = zeros(1,200);
temp_linear = sum(genome_gen_linear{gen},2);
temp_linear = reshape(temp_linear, [1000 200]);
is_disease_resist_linear = find(temp_linear(118,:) <= 1 & temp_linear(533,:) <= 1 | temp_linear(831,:) > 1);
disease_pheno_linear(1, is_disease_resist_linear) = 1;


disease_pheno_nonlinear = zeros(1,200);
temp_nonlinear = sum(genome_gen_nonlinear{gen},2);
temp_nonlinear = reshape(temp_nonlinear, [1000 200]);
is_disease_resist_nonlinear = find(temp_nonlinear(118,:) <= 1 & temp_nonlinear(533,:) <=1 | temp_nonlinear(831,:) > 1);
disease_pheno_nonlinear(1, is_disease_resist_nonlinear) = 1;

end


function objective_function_nonlinear = f_obj_nonlinear(pheno_200,minh, maxh, rvalue,gen)
% N = 10;


%% normalizing traits

pheno_200(1,:) = (pheno_200(1,:)-minh)/(maxh-minh);
% pheno_200(2,:) = (pheno_200(2,:)-minh)/(maxh-minh);
pheno_200(2,:) = (pheno_200(2,:)) ;
pheno_200(pheno_200 == 0 ) = 0.01; 
% pheno_200(2,:) = rescale(pheno_200(2,:), 0.01 ,1);

a = find(pheno_200 > 1 | pheno_200 <= 0 );
if ( isempty(a) ~= 1)
    disp('nonlinear out of bound');
end
% pheno_180(3,:) = pheno_180(3,:)/20;

temp = [(pheno_200(1,:)/rvalue(1)) ; (pheno_200(2,:)/rvalue(2))];
%% calculating objective function
objective_function_nonlinear = min(temp);
end


function objective_function_linear = f_obj_linear(pheno_200,minh, maxh, rvalue,gen)


pheno_200(1,:) = (pheno_200(1,:)-minh)/(maxh-minh);
% pheno_200(2,:) = (pheno_200(2,:)-minh)/(maxh-minh);
pheno_200(2,:) = (pheno_200(2,:));

pheno_200(pheno_200 == 0 ) = 0.01;
% pheno_200(2,:) = rescale(pheno_200(2,:), 0.01, 1) ;
% pheno_180(3,:) = pheno_180(3,:)/20;
a = find(pheno_200 > 1 | pheno_200 <= 0  );
if ( isempty(a) ~= 1)
    disp('linear out of bound');
end
%% calculating objective function
objective_function_linear =  rvalue(1).*pheno_200(1,:) +...
    rvalue(2).*pheno_200(2,:);
end




% select 20 best parents in each generation.

function final_parents_nonlinear = find_parents_nonlinear(pheno_all_nonlinear, gen, minh, maxh, r, ri)
    
    pheno_nonlinear = pheno_all_nonlinear(2:end, :);
    pheno_200 = f_obj_nonlinear(pheno_nonlinear, minh, maxh, r{ri},gen);
     [~, final_parents_nonlinear] = maxk(pheno_200, 40);
    
end

% same function as the above function for linear algorithm
function final_parents_linear = find_parents_linear( pheno_all_linear, gen , minh, maxh, r , ri)

    pheno_linear = pheno_all_linear(2:end, :);
    pheno_200 = f_obj_linear(pheno_linear, minh, maxh,  r{ri},gen);
     [~, final_parents_linear] = maxk(pheno_200, 40);

end



function [pop_linear, pop_nonlinear] = cross(genome_gen_linear,genome_gen_nonlinear, final_parents_linear, final_parents_nonlinear, recomData,   gen)

pool_cross_linear = randperm(40);
pool_cross_nonlinear = randperm(40);


pop_linear =[];
pop_nonlinear = [];
j = 1;
while (j < 20) % numebr of crosses = 10/2 = 5
    %
    k = 20; % number of progeny per cross
    L1_linear = genome_gen_linear{gen}(:,:,final_parents_linear(pool_cross_linear(j)));
    L2_linear = genome_gen_linear{gen}(:,:,final_parents_linear(pool_cross_linear(j+1)));
    L1_nonlinear = genome_gen_nonlinear{gen}(:,:,final_parents_nonlinear(pool_cross_nonlinear(j)));
    L2_nonlinear = genome_gen_nonlinear{gen}(:,:,final_parents_nonlinear(pool_cross_nonlinear(j+1)));
    
    n1 = size(L1_linear, 1);
    n2 = size(L1_nonlinear, 1);
    
    RC_linear = rand(n1, 2*k) <= repmat([0.5; recomData], 1, 2*k);
    RC_nonlinear = rand(n2, 2*k) <= repmat([0.5; recomData], 1, 2*k);
    
    
    
    
    Y1_linear = repmat([L1_linear(:, 1), L2_linear(:, 1)], 1, k);
    Y2_linear = repmat([L1_linear(:, 2), L2_linear(:, 2)], 1, k);
    Y1_nonlinear = repmat([L1_nonlinear(:, 1), L2_nonlinear(:, 1)], 1, k);
    Y2_nonlinear = repmat([L1_nonlinear(:, 2), L2_nonlinear(:, 2)], 1, k);
    
    f_linear = cumprod(1 - 2 * RC_linear) <= 0;
    f_nonlinear = cumprod(1 - 2 * RC_nonlinear) <= 0;
    
    Y1_linear(f_linear) = Y2_linear(f_linear);
    Y1_nonlinear(f_nonlinear) = Y2_nonlinear(f_nonlinear);
    
    Y1_linear = reshape(Y1_linear, n1, 2, k);
    Y1_nonlinear = reshape(Y1_nonlinear, n2, 2, k);
    
    pop_linear = cat(3,pop_linear,Y1_linear);
    pop_nonlinear = cat(3,pop_nonlinear, Y1_nonlinear);
    
    j = j + 2;
    
end
end



% this function is used to create \beta values

function [x,v] = randfixedsum(n,m,s,a,b)

% [x,v] = randfixedsum(n,m,s,a,b)
%
%   This generates an n by m array x, each of whose m columns
% contains n random values lying in the interval [a,b], but
% subject to the condition that their sum be equal to s.  The
% scalar value s must accordingly satisfy n*a <= s <= n*b.  The
% distribution of values is uniform in the sense that it has the
% conditional probability distribution of a uniform distribution
% over the whole n-cube, given that the sum of the x's is s.
%
%   The scalar v, if requested, returns with the total
% n-1 dimensional volume (content) of the subset satisfying
% this condition.  Consequently if v, considered as a function
% of s and divided by sqrt(n), is integrated with respect to s
% from s = a to s = b, the result would necessarily be the
% n-dimensional volume of the whole cube, namely (b-a)^n.
%
%   This algorithm does no "rejecting" on the sets of x's it
% obtains.  It is designed to generate only those that satisfy all
% the above conditions and to do so with a uniform distribution.
% It accomplishes this by decomposing the space of all possible x
% sets (columns) into n-1 dimensional simplexes.  (Line segments,
% triangles, and tetrahedra, are one-, two-, and three-dimensional
% examples of simplexes, respectively.)  It makes use of three
% different sets of 'rand' variables, one to locate values
% uniformly within each type of simplex, another to randomly
% select representatives of each different type of simplex in
% proportion to their volume, and a third to perform random
% permutations to provide an even distribution of simplex choices
% among like types.  For example, with n equal to 3 and s set at,
% say, 40% of the way from a towards b, there will be 2 different
% types of simplex, in this case triangles, each with its own
% area, and 6 different versions of each from permutations, for
% a total of 12 triangles, and these all fit together to form a
% particular planar non-regular hexagon in 3 dimensions, with v
% returned set equal to the hexagon's area.
%
% Roger Stafford - Jan. 19, 2006

% Check the arguments.
if (m~=round(m))|(n~=round(n))|(m<0)|(n<1)
    error('n must be a whole number and m a non-negative integer.')
elseif (s<n*a)|(s>n*b)|(a>=b)
    error('Inequalities n*a <= s <= n*b and a < b must hold.')
end

% Rescale to a unit cube: 0 <= x(i) <= 1
s = (s-n*a)/(b-a);

% Construct the transition probability table, t.
% t(i,j) will be utilized only in the region where j <= i + 1.
k = max(min(floor(s),n-1),0); % Must have 0 <= k <= n-1
s = max(min(s,k+1),k); % Must have k <= s <= k+1
s1 = s - [k:-1:k-n+1]; % s1 & s2 will never be negative
s2 = [k+n:-1:k+1] - s;
w = zeros(n,n+1); w(1,2) = realmax; % Scale for full 'double' range
t = zeros(n-1,n);
tiny = 2^(-1074); % The smallest positive matlab 'double' no.
for i = 2:n
    tmp1 = w(i-1,2:i+1).*s1(1:i)/i;
    tmp2 = w(i-1,1:i).*s2(n-i+1:n)/i;
    w(i,2:i+1) = tmp1 + tmp2;
    tmp3 = w(i,2:i+1) + tiny; % In case tmp1 & tmp2 are both 0,
    tmp4 = (s2(n-i+1:n) > s1(1:i)); % then t is 0 on left & 1 on right
    t(i-1,1:i) = (tmp2./tmp3).*tmp4 + (1-tmp1./tmp3).*(~tmp4);
end

% Derive the polytope volume v from the appropriate
% element in the bottom row of w.
v = n^(3/2)*(w(n,k+2)/realmax)*(b-a)^(n-1);

% Now compute the matrix x.
x = zeros(n,m);
if m == 0, return, end % If m is zero, quit with x = []
rt = rand(n-1,m); % For random selection of simplex type
rs = rand(n-1,m); % For random location within a simplex
s = repmat(s,1,m);
j = repmat(k+1,1,m); % For indexing in the t table
sm = zeros(1,m); pr = ones(1,m); % Start with sum zero & product 1
for i = n-1:-1:1  % Work backwards in the t table
    e = (rt(n-i,:)<=t(i,j)); % Use rt to choose a transition
    sx = rs(n-i,:).^(1/i); % Use rs to compute next simplex coord.
    sm = sm + (1-sx).*pr.*s/(i+1); % Update sum
    pr = sx.*pr; % Update product
    x(n-i,:) = sm + pr.*e; % Calculate x using simplex coords.
    s = s - e; j = j - e; % Transition adjustment
end
x(n,:) = sm + pr.*s; % Compute the last x

% Randomly permute the order in the columns of x and rescale.
rp = rand(n,m); % Use rp to carry out a matrix 'randperm'
[ig,p] = sort(rp); % The values placed in ig are ignored
x = (b-a)*x(p+repmat([0:n:n*(m-1)],n,1))+a; % Permute & rescale x

return
end




function [minor_allele_linear, minor_allele_nonlinear] = minor_allele_calc(genome_gen_linear, genome_gen_nonlinear, gen)
%calculate minor alles

%     if (gen == 1)
zeros_sum_linear = nnz(~genome_gen_linear{gen})/(1000*2*200);
ones_sum_linear = nnz(genome_gen_linear{gen})/(1000*2*200) ;

zeros_sum_nonlinear = nnz(~genome_gen_nonlinear{gen}) /(1000*2*200) ;
ones_sum_nonlinear = nnz(genome_gen_nonlinear{gen})/(1000*2*200) ;
%     else
%         zeros_sum_linear = nnz(~genome_gen_linear{gen})/(1000*2*1000);
%         ones_sum_linear = nnz(genome_gen_linear{gen})/(1000*2*1000) ;
%
%         zeros_sum_nonlinear = nnz(~genome_gen_nonlinear{gen}) /(1000*2*1000) ;
%         ones_sum_nonlinear = nnz(genome_gen_nonlinear{gen})/(1000*2*1000) ;
%
%     end

if (zeros_sum_linear > ones_sum_linear)
    minor_allele_linear = ones_sum_linear;
else
    minor_allele_linear = zeros_sum_linear;
end

%     minor_allele_linear_all = [minor_allele_linear_all;minor_allele_linear ];


if (zeros_sum_nonlinear > ones_sum_nonlinear)
    minor_allele_nonlinear = ones_sum_nonlinear;
else
    minor_allele_nonlinear = zeros_sum_nonlinear;
end

%     minor_allele_nonlinear_all = [minor_allele_nonlinear_all;minor_allele_nonlinear ];
end

