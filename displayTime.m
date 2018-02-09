function displayTime(secs)

fprintf('Total time = %5.0f seconds, or %5.0f minutes and %5.0f seconds\n',...
    secs, floor(secs/60), secs - floor(secs/60) * 60)