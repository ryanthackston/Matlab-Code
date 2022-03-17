function [MIcon_Imag_Aug2021] = micon_imag_aug2021()
%Initial inputs for MEG study
%   Detailed explanation goes here
% Load data and rearrange the values
% Cell 1 is Subj 12
% Cells 2-5 is Subj 13
MIcon_Imag_Aug2021 ={ { { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block4_imag_control.mat')] } };...
                                      
                                      
                                   {  { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D2_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D2_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D2_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D2_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D4_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D4_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D4_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D4_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D5_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D5_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D5_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S13_D5_Block4_imag_control.mat')] } };...
                                      
                                      
                                  {   { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D2_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D2_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D2_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D2_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D3_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D3_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D3_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D3_Block4_imag_control.mat')] };
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D4_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D4_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D4_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D4_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D5_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D5_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S14_D5_Block3_imag_control.mat')] } };...
                                      
                                      
                                  {   { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_Block3_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D2_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D2_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D2_Block3_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D3_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D3_Block2_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D3_Block4_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D4_Block1_imag_control.mat')] };...
                                      
                                      { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S15_D5_Block1_imag_control.mat')] } };...
                                      
                                      
                                 {    { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S16_Block1_imag_control.mat');...
                                      ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S16_Block4_imag_control.mat')] } };...
                                      
                                      
                                 {    { [('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S18_Block2_imag_control.mat')] } } };
                                  
                                  
%                                       {[('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block1_imag_control.mat');...
%                                       ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block2_imag_control.mat');...
%                                       ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block3_imag_control.mat');...
%                                       ('C:\Users\ryanr\Downloads\My Matlab Functions\OPM data - August 2021\Control MAT files\S12_Block4_imag_control.mat')]}

%                                     subj = 13;
%                                     z = [];
%                                     for D = 1: size(MIcal_Imag_Aug2021{subj-11}{1}{1},1)
%                                         tmp = MIcal_Imag_Aug2021{subj-11}{D}{1};
%                                         z = char(z,tmp);
%                                     end
%                                     z(1,:) = [];
                                       

end
