function r=integrand_2_on_edges(x,y,basis_type,derivative_degree_x,derivative_degree_y)
%The integrand of the integral on the inner edges

r=triangular_local_basis(x,y,basis_type,derivative_degree_x,derivative_degree_y)'...
    *triangular_local_basis(x,y,basis_type,derivative_degree_x,derivative_degree_y);

