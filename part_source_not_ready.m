function [ part_source ] = part_source( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



part_source(:,:,1) = [tvect' [zeros(1,66), 1.0.*ones(1,12), zeros(1,67)]' (3e-9.*ones(1,length(tvect)))'];

end

