%% Renyi dimension of Hénon map.
% Parameters.
clearvars
tic
a = 1.4; b = 0.3; nMax = 800; sqrtInit = 50; numberOfBoxes = 10^3;

% Initialization.
x0 = linspace(-0.1, 0.1, sqrtInit); y0 = linspace(-0.1, 0.1, sqrtInit);
[x0s, y0s] = meshgrid(x0, y0);
xList = zeros(sqrtInit, sqrtInit, nMax);
yList = zeros(sqrtInit, sqrtInit, nMax);

% Main loop.
for i = 1:sqrtInit
    for j = 1:sqrtInit
        x = zeros(1, nMax);
        y = zeros(1, nMax);
        x(1) = x0s(i, j);
        y(1) = y0s(i, j);
        for n = 1:(nMax - 1)
            x(n + 1) = y(n) + 1 - a*x(n).^2;
            y(n + 1) = b*x(n);
        end
        xList(i, j, :) = x(:);
        yList(i, j, :) = y(:);
    end
end

% Make boxes.
xMin = min(min(min(xList))); xMax = max(max(max(xList)));
%yMin = min(min(min(yList))); yMax = max(max(max(yList)));
boxXs = linspace(xMin, xMax, numberOfBoxes); % Wanted boxes to be square.
boxYs = linspace(xMin, xMax, numberOfBoxes);
boxesCounter = zeros(numberOfBoxes, numberOfBoxes);
for i = 1:sqrtInit
    for j = 1:sqrtInit
        for n = 25:nMax
            x = xList(i, j, n);
            y = yList(i, j, n);
            for m = 2:numberOfBoxes
                if (~(isempty(find(x >= boxXs(m - 1) & x <= boxXs(m), 1))))
                    for k = 2:numberOfBoxes
                        if (~isempty(find(y >= boxYs(k - 1) & y <= boxYs(k), 1)))
                            boxesCounter(k, m) = boxesCounter(k, m) + 1;
                            break
                        end
                    end
                end
            end
        end
    end
    disp(i)
end
save('boxesCounter', 'boxesCounter');

%% Calculation of Dq.
Q = [0, 2];
Iq = [0, 0];
Ntot = nMax*sqrtInit^2;
for q = 1:length(Q)
    for j = 1:sqrtInit^2
        Iq(q) = Iq(q) + (boxesCounter(j)/Ntot)^Q(q);
    end
end

epsilon = (xMax - xMin) / numberOfBoxes;
Dq = zeros(1, 3);
for q = 1:length(Q)
    Dq(q) = 1/(1 - Q(q))*log(Iq(q)) / log(1/epsilon);
end

for k = 1:sqrtInit^2
    Dq(end) = Dq(end) + (boxesCounter(j)/Ntot)*log(boxesCounter(j)/Ntot)/log(epsilon);
end

toc

%% Plot of the points.
hold on
title('Approximation of the Hénon map.')
for i = 1:sqrtInit
    for j = 1:sqrtInit
        plot(squeeze(xList(i, j, 25:end)), squeeze(yList(i, j, 25:end)), '.')
    end
end
hold off

%% Plot of the boxes.
load('boxesCounter', 'boxesCounter')
set(gca, 'YDir', 'normal')
colormap(flipud(gray))
imagesc(boxesCounter)
colorbar
pbaspect([1 1 1])
xlabel('x')
ylabel('y')
title("Plot of the boxes.")
