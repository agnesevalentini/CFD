% definition of u for Poiseuille

function my_u = u(y, nu, b, dpdx)
  my_u = -((y ./ nu) .* dpdx) .* (b - (y ./ 2));
 endfunction
