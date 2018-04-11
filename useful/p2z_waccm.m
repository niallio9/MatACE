function [ z ] = p2z_waccm( p )
%pressure is in Pa and z is in m

z = -7000*log(p./1e5);

end

