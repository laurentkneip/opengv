function R = cayley2rot(v)

  cayley0 = v(1,1);
  cayley1 = v(2,1);
  cayley2 = v(3,1);
  
  R = zeros(3,3);
  scale = 1+cayley0^2+cayley1^2+cayley2^2;

  R(1,1) = 1+cayley0^2-cayley1^2-cayley2^2;
  R(1,2) = 2*(cayley0*cayley1-cayley2);
  R(1,3) = 2*(cayley0*cayley2+cayley1);
  R(2,1) = 2*(cayley0*cayley1+cayley2);
  R(2,2) = 1-cayley0^2+cayley1^2-cayley2^2;
  R(2,3) = 2*(cayley1*cayley2-cayley0);
  R(3,1) = 2*(cayley0*cayley2-cayley1);
  R(3,2) = 2*(cayley1*cayley2+cayley0);
  R(3,3) = 1-cayley0^2-cayley1^2+cayley2^2;

  R = (1/scale) * R;