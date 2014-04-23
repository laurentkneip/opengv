function v = rot2cayley(R)

  C1 = R-eye(3);
  C2 = R+eye(3);
  C = C1 * inv(C2);

  v = zeros(3,1);
  v(1,1) = -C(2,3);
  v(2,1) =  C(1,3);
  v(3,1) = -C(1,2);

end