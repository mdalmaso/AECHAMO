function[out_struct] = model_convert(obj,t,Y)

if(obj.initials.fixed_sections == 0)
    out_struct = model_convert_moving_sec(obj,t,Y);
else
%     out_struct.orig = model_convert_moving_sec(obj,t,Y);
%     out_struct.new = model_convert_moving_center(obj,t,Y);
    out_struct = model_convert_moving_center(obj,t,Y);
end

end