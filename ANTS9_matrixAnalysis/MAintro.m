%% Mike X Cohen

%% vectors as lines

% 2-dimensional vector
v2 = [ 3 -2 ];

% 3-dimensional vector
v3 = [ 4 -3 2 ];


% plot them
figure(1), clf
subplot(211)
plot([0 v2(1)],[0 v2(2)],'linew',2)
axis square
axis([ -4 4 -4 4 ])
hold on
plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')


subplot(212)
plot3([0 v3(1)],[0 v3(2)],[0 v3(3)],'linew',2)
axis square
axis([ -4 4 -4 4 -4 4 ])
hold on, grid on
plot3(get(gca,'xlim'),[0 0],[0 0],'k--')
plot3([0 0],get(gca,'ylim'),[0 0],'k--')
plot3([0 0],[0 0],get(gca,'zlim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')
zlabel('X_3 dimension')

% might be easier to see when rotated
rotate3d on

%% vector-vector addition

v1 = [3 2];
v2 = [1 2];
v3 = v1+v2;

figure(2), clf
subplot(121)
plot([0 v1(1)],[0 v1(2)],'k','linew',2)
axis square
axis([ -6 6 -6 6 ])
hold on
plot([0 v2(1)],[0 v2(2)],'b','linew',2)
plot([0 v3(1)],[0 v3(2)],'r','linew',2)

plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')

legend({'v1';'v2';'v1+v2'},'location','northwest') % also specify the location of the legend

%% vector-scalar multiplication

v2m = v2 * 1.5;

% plot them
subplot(122)
plot([0 v2(1)],[0 v2(2)],'k','linew',2)
axis square
axis([ -6 6 -6 6 ])
hold on
plot([0 v2m(1)],[0 v2m(2)],'r--','linew',2)

plot(get(gca,'xlim'),[0 0],'k--')
plot([0 0],get(gca,'ylim'),'k--')
xlabel('X_1 dimension')
ylabel('X_2 dimension')

legend({'v2';'1.5v2'})

%% matrix addition and scalar multiplication

% matrix addition requires two matrices of the same size
A1 = round( 10*randn(2,3) );
A2 = round( 10*randn(2,3) );
A3 = round( 10*randn(3,2) );

A1 + A2
A2 + A1 % addition is commutative

A1 + A3


% scalar multiplication always works
4*A1
.01*A2
A3*10

% of course, you can combine scalar multiplication and addition
C = 2.5*A1 - 6.2*A2

%% rules for matrix multiplication

A = randn(2,4);
B = randn(3,2);
C = randn(6,4);

A*B
B*A

%% finding your way around matrices

% create a matrix by reshaping a vector
A = reshape( 1:12, 3,4)

% matrix indexing: state the element you want in the dimensions you want
A(3,2)
A(2:3,3:4)
A([1 3],[1 3 4])

% linear indexing: use one number to access elements (see figure 10.6)
A(4)
A(11)
A(3:5)
A([2 4 10 12])


% the same principle applies in higher dimensions
B = randn(4,2,6,5,3);

% please please please do yourself a favor and 
% use matrix indexing whenever possible!
B(1,2,3,2,1)

% linear indexing can get quite confusing. 
% Where in the matrix is this element??
B(582)


%%

%% transposing vectors

v = rand(10,1); % 10-dimensional column vector
v
size(v)
v'
size(v')
transpose(v)

% transpose of a transpose is the original vector
v==v''
% question: why does the previous line return 10 "1"s?

%% Three important matrices

N = 4; % for easy viewing and inspection

% square matrix
sqr = round( 10*rand(N) );

% symmetric matrix
% Use the rule that a matrix times its transpose is symmetric.
sym = round( sqr'*sqr ); 

% identity matrix
idt = eye(N);

figure(3), clf
subplot(131), imagesc(sqr), axis square, title('Square')
subplot(132), imagesc(sym), axis square, title('Symmetric')
subplot(133), imagesc(idt), axis square, title('Identity')

%% rank

% random matrices are basically always full-rank
rank( rand(4)  )
rank( rand(4,3) )

% now for a reduced-rank matrix
A = [1 2 3; 2 4 6; 3 6 8];
rank(A)
% A is a 3x3 matrix, but has rank=2.

%% end.
