% Make a movie of N convolutions and associated gaussians for a few
% different function distributions.
% 
% INPUT:
%           N: number of convolutions
% 
%
%  functions used:
%       exponential
%       polynomial
%       sin (made positive)
%
%
%   last edited by lrgrant 4 October 2022

N = 24;
convolution_movie(N)

function convolution_movie(N)
% an x space
x = linspace(-4, 4, 100);

% a few random functions, made into distribution functions [a, b, c]
a = exp(x);
a_n = trapz(x,a);
a = a/a_n;

b = sin(x) + 1;
b_n = trapz(x,b);
b = b/b_n;

c = 7*x.^2 + 5*x - 8;
c_n = trapz(x,c);
c = c/c_n;

% perform first convolution
[x, y] = convolve(a,b);

% find mean and var
mean = trapz(x,x.*y);
var = trapz(x,(x-mean).^2.*y);

% make gaussian of same mean and var
normal = normpdf(x, mean, sqrt(var));
xlim(mean+[-5 5]*sqrt(var));

% save each convolution and associated gaussian as a movie frame
% create movie
M(N) = struct('cdata',[],'colormap',[]);

% create and hide figure
f = figure;
f.Visible = 'off';

% plot and save convolutions
plot(x, y, linewidth = 2)
hold on
plot(x, normal, linewidth = 2)
legend("Convolved", "Gaussian", location='northwest')
set(gca, fontsize=16)
M(1) = getframe;

for i = 2:N-1
    % convolve with either a, b, or c
    if mod(i,3) == 1
        [x,y] = convolve(a,y);
    elseif mod(i,3) == 2
        [x,y] = convolve(b,y);
    else
        [x,y] = convolve(c,y);
    end

    % find mean and var
    mean = trapz(x,x.*y);
    var = trapz(x,(x-mean).^2.*y);

    % make gaussian of same mean and var
    normal = normpdf(x, mean, sqrt(var));
    xlim(mean+[-5 5]*sqrt(var));

    % plot and save convolution and gaussian as movie frame
    clf;
    plot(x, y, LineWidth=2)
    hold on
    plot(x, normal, linewidth = 2)
    legend("Convolved", "Gaussian", location='northwest')
    set(gca, fontsize=16)
    M(i) = getframe;
    clf;
end

% play movie
figure(1)
ylabel("Frequency")
xlabel("X")
title(sprintf("%i Successive Convolutions", N))
set(gca, fontsize=16)
movie(M, 1, 5)


% my function to convolve and normalize
function [x,y] = convolve(a,b)
% convolve
y = conv(a,b);
% readjust x space
x = linspace(-4,4,length(y));
% normalize
n = trapz(x,y);
y = y/n;
end
end
