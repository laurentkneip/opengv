function R = generateBoundedR( bound )
  rpy = bound*2.0*(rand(3,1)-repmat(0.5,3,1));

  R1 = zeros(3,3);
  R1(1,1) = 1.0;
  R1(2,2) = cos(rpy(1,1));
  R1(2,3) = -sin(rpy(1,1));
  R1(3,2) = -R1(2,3);
  R1(3,3) = R1(2,2);

  R2 = zeros(3,3);
  R2(1,1) = cos(rpy(2,1));
  R2(1,3) = sin(rpy(2,1));
  R2(2,2) = 1.0;
  R2(3,1) = -R2(1,3);
  R2(3,3) = R2(1,1);

  R3 = zeros(3,3);
  R3(1,1) = cos(rpy(3,1));
  R3(1,2) = -sin(rpy(3,1));
  R3(2,1) =-R3(1,2);
  R3(2,2) = R3(1,1);
  R3(3,3) = 1.0;

  R = R3 * R2 * R1;
end