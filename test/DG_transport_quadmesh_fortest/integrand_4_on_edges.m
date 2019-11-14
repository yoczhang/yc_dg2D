function r=integrand_4_on_edges(x,y,basis_type)

r=quadrilateral_local_basis(x,y,basis_type,0,0)'...
    *quadrilateral_local_basis(x,y,basis_type,0,0);