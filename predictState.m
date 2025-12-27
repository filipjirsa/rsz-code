function [x_pred, F] = predictState(x, acc, gyr, dt)

p  = x(1:3); v  = x(4:6);
q  = x(7:10);
ba = x(11:13); bg = x(14:16);

acc = acc - ba;
gyr = gyr - bg;

R = quat2rotm(q');
a_world = R * acc;

% stav
p = p + v*dt + 0.5*a_world*dt^2;
v = v + a_world*dt;

% quaternion
wx=gyr(1); wy=gyr(2); wz=gyr(3);
Om = [0 -wx -wy -wz;
      wx 0 wz -wy;
      wy -wz 0 wx;
      wz wy -wx 0];

q = (eye(4)+0.5*Om*dt)*q;
q = q / norm(q);

x_pred = [p; v; q; ba; bg];

% Jacobián (zjednodušený, ale konzistentný)
F = eye(16);
F(1:3,4:6) = eye(3)*dt;
F(4:6,11:13) = -R*dt;
F(7:10,7:10) = eye(4) + 0.5*Om*dt;
end
