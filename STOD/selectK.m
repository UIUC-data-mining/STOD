function K = selectK(eigenvals, options)
upper = length(eigenvals);
lower = min(options.K);
if upper <= lower
  K = upper;
  return;
end
sumofall = sum(eigenvals);
for K = options.K
    sumoffirstK = sum(eigenvals(1:K));
    fprintf('first %d: %f of %d\n', K, sumoffirstK/sumofall, K);
    if sumoffirstK >= options.proportion * sumofall
      break;
    end
end  
