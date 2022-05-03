a= T.Duration;
b= T.LPON_MeanAbs;
% First need to sort a otherwise when we go to plot it, it will look like a mess!
[~, sortOrder] = sort(a, 'ascend');
a = a(sortOrder);
b = b(sortOrder); % Need to sort b the same way.
indices = find(isnan(b) == 1);
a(indices) = [];
b(indices) = [];
% First compute the linear fit.
linearCoeffs = polyfit(a, b, 1);
Slope = linearCoeffs(2)
Intercept = linearCoeffs(1)
% Plot training data and fitted data.
subplot(2, 1, 1);
aFitted = a; % Evalutate the fit as the same x coordinates.
bFitted = polyval(linearCoeffs, aFitted);
plot(a, b, 'rd', 'MarkerSize', 10);
hold on;
plot(aFitted, bFitted, 'b-', 'LineWidth', 2);
grid on;
xlabel('a', 'FontSize', 20);
ylabel('b', 'FontSize', 20);
% Plot residuals as lines from actual data to fitted line.
for k = 1 : length(a)
  yActual = b(k);
  yFit = bFitted(k);
  x = a(k);
  line([x, x], [yFit, yActual], 'Color', 'm');
end
% Do the same for a quadratic fit.
quadraticCoeffs = polyfit(a, b, 2);
% Plot training data and fitted data.
subplot(2, 1, 2);
aFitted = a; % Evalutate the fit as the same x coordinates.
bFitted = polyval(quadraticCoeffs, aFitted);
plot(a, b, 'rd', 'MarkerSize', 10);
hold on;
plot(aFitted, bFitted, 'b-', 'LineWidth', 2);
grid on;
xlabel('a', 'FontSize', 20);
ylabel('b', 'FontSize', 20);
% Plot residuals as lines from actual data to fitted line.
for k = 1 : length(a)
  yActual = b(k);
  yFit = bFitted(k);
  x = a(k);
  line([x, x], [yFit, yActual], 'Color', 'm');
end