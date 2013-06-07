function run(obj)
if(obj.initials.fixed_sections == 0)
%     tic
    obj.run_movsec;
%     toc
else
%     tic
    obj.run_moving_center;
%     toc
end

end