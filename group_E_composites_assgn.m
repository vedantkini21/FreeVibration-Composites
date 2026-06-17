%Levy solution for SSSF boundary conditions
%Laminate: 90/0/90 degrees

clear all; close all; clc;

%% material properties - given in the assignment
E11 = 18;%E1 is 18 times E2
E22 = 1;%using E2 as reference value
G12 = 0.8;%shear modulus
v12 = 0.28;%Poisson's ratio
%calculate the other Poisson's ratio
v21 = v12 * (E22/E11);

%calculate the stiffness matrix components
Q11 = E11/(1 - v12*v21);
Q22 = E22/(1 - v12*v21);
Q12 = v12*E22/(1 - v12*v21);
Q66 = G12;

%build the full Q matrix
Q = [Q11, Q12, 0;
     Q12, Q22, 0;
     0,   0, Q66];

%% laminate setup - three layers 90/0/90
n_layers = 3;
h_total = 1;%total thickness
h_layer = h_total/n_layers;%each layer same thickness

%layer angles in degrees
theta = [90, 0, 90];

%initialize the ABD matrices
A = zeros(3); B = zeros(3); D = zeros(3);
z0 = -h_total/2;  % bottom of laminate

%loop through each layer to calculate ABD matrices
for k = 1:n_layers
    z_k_minus = z0 + (k-1)*h_layer;
    z_k = z0 + k*h_layer;
    
    %transformation matrix for this layer's angle
    angle_rad = theta(k) * pi/180;
    m = cos(angle_rad);
    n = sin(angle_rad);
    
    T = [m^2, n^2, 2*m*n;
         n^2, m^2, -2*m*n;
         -m*n, m*n, m^2-n^2];
    
    %transformed stiffness matrix
    Q_bar = T \ Q / T.';
    
    %ABD matrices
    A = A + Q_bar * (z_k - z_k_minus);
    B = B + 0.5 * Q_bar * (z_k^2 - z_k_minus^2);
    D = D + (1/3) * Q_bar * (z_k^3 - z_k_minus^3);
end

%extracting the bending stiffness values
D11 = D(1,1); D12 = D(1,2); D22 = D(2,2); D66 = D(3,3);
D16 = D(1,3); D26 = D(2,3);

fprintf('bending stiffness values:\n');
fprintf('D11 = %.4f, D12 = %.4f, D22 = %.4f\n', D11, D12, D22);
fprintf('D66 = %.4f\n\n', D66);

%% problem 1 - first two frequencies for square plate
a = 1;%plate length
b = 1; %plate width  
rho = 1;%density

%for SSSF simply supported at x=0,a and y=0; free at y=b

max_m = 8;%number of terms in x-direction
max_n = 6;%number of terms for frequency calculation

%storing frequencies
all_frequencies = [];

%looping through different mode numbers
for m = 1:max_m
    alpha = m * pi / a;
    
    %for each m, check different n values
    for n = 1:max_n
        %approximate wave numbers for free edge
        if n == 1
            beta = 1.875 / b;  % first mode for free edge
        elseif n == 2
            beta = 4.694 / b;  % second mode
        else
            beta = ((2*n-1) * pi) / (2 * b);  % higher modes
        end
        
        %calculating frequency parameter
        lambda = D11*alpha^4 + 2*(D12 + 2*D66)*alpha^2*beta^2 + D22*beta^4;
        omega = sqrt(lambda / rho);
        
        %storing frequency with mode numbers
        all_frequencies = [all_frequencies; omega, m, n];
    end
end

%sort by frequency to find lowest ones
all_frequencies = sortrows(all_frequencies, 1);
omega1 = all_frequencies(1,1);
omega2 = all_frequencies(2,1);
m1 = all_frequencies(1,2); n1 = all_frequencies(1,3);
m2 = all_frequencies(2,2); n2 = all_frequencies(2,3);

fprintf('first frequency: %.4f rad/s (m=%d, n=%d)\n', omega1, m1, n1);
fprintf('second frequency: %.4f rad/s (m=%d, n=%d)\n\n', omega2, m2, n2);

%% plots of the mode shapes - first two modes

%create grid for plotting
nx = 50; ny = 50;
x = linspace(0, a, nx);
y = linspace(0, b, ny);
[X, Y] = meshgrid(x, y);

%initialize mode shape arrays
W1 = zeros(ny, nx);
W2 = zeros(ny, nx);

%calculate mode shapes point by point
for i = 1:nx
    for j = 1:ny
        %x-direction is sine function for simply supported
        X_shape1 = sin(m1 * pi * x(i) / a);
        X_shape2 = sin(m2 * pi * x(i) / a);
        
        %y-direction shape - approximate for free edge
        Y_shape1 = sin(pi * y(j) / (2*b)) + 0.5 * sin(3*pi * y(j) / (2*b));
        Y_shape2 = sin(3*pi * y(j) / (2*b)) + 0.3 * sin(5*pi * y(j) / (2*b));
        
        W1(j,i) = X_shape1 * Y_shape1;
        W2(j,i) = X_shape2 * Y_shape2;
    end
end

%normalize so max is 1
W1 = W1 / max(abs(W1(:)));
W2 = W2 / max(abs(W2(:)));

%plots of both mode shapes
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1);
surf(X, Y, W1, 'EdgeColor', 'none');
title(sprintf('First Mode Shape (ω=%.4f)', omega1));
xlabel('x'); ylabel('y'); zlabel('w');
colorbar; axis equal;

subplot(1,2,2);
surf(X, Y, W2, 'EdgeColor', 'none');
title(sprintf('Second Mode Shape (ω=%.4f)', omega2));
xlabel('x'); ylabel('y'); zlabel('w');
colorbar; axis equal;

%% problem 2 - frequency vs aspect ratio
aspect_ratios = 1:0.2:5;%from 1 to 5
first_frequencies = zeros(size(aspect_ratios));

%lowest frequency for each aspect ratio
for idx = 1:length(aspect_ratios)
    ar = aspect_ratios(idx);
    a_current = ar;%change of length
    b_current = 1;%width kept fixed
    
    %find minimum frequency for this aspect ratio
    omega_min = inf;
    best_m = 1;
    
    for m = 1:6
        alpha = m * pi / a_current;
        beta = 1.875 / b_current; %fundamental free edge mode
        
        lambda = D11*alpha^4 + 2*(D12 + 2*D66)*alpha^2*beta^2 + D22*beta^4;
        omega_current = sqrt(lambda / rho);
        
        if omega_current < omega_min
            omega_min = omega_current;
            best_m = m;
        end
    end
    
    first_frequencies(idx) = omega_min;
    fprintf('a/b=%.1f: frequency=%.4f (m=%d)\n', ar, omega_min, best_m);
end

%normalizing by the value at a/b=1
freq_normalized = first_frequencies / first_frequencies(1);

%plots
figure;
plot(aspect_ratios, freq_normalized, 'r-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
xlabel('aspect ratio (a/b)');
ylabel('normalized frequency');
title('First Frequency vs Aspect Ratio');
grid on;

%display results table
fprintf('\nsummary of aspect ratio study:\n');
fprintf('a/b\tnormalized frequency\n');
fprintf('----------------------------\n');
for i = 1:length(aspect_ratios)
    fprintf('%.1f\t%.4f\n', aspect_ratios(i), freq_normalized(i));
end

%% final summary
fprintf('laminate: 90/0/90 degrees\n');
fprintf('boundary conditions: SSSF\n');
fprintf('square plate results:\n');
fprintf('  first frequency:  %.4f rad/s\n', omega1);
fprintf('  second frequency: %.4f rad/s\n', omega2);



%ABD matrices
fprintf('\nABD matrices:\n');
fprintf('A matrix:\n');
disp(A);
fprintf('B matrix:\n');
disp(B);
fprintf('D matrix:\n');
disp(D);