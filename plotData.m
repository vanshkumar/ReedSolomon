function plotData(name)

fid = fopen(name);
s = textscan(fid, '%f %f %f');
probs = fliplr(s{1});
ser = fliplr(s{2});
wer = fliplr(s{3});

plot(probs, ser, 'Color', 'red');
hold on;
plot(probs, wer, 'Color', 'blue')
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
legend('Symbol Error Rate', 'Word Error Rate')

end