function [high_pass_signal, filtered_signal] = high_pass_filter_data(signal, std_dev)


%high_pass_signal = zeros(size(signal));



%std_dev = 20;
if ~exist('std_dev')
    std_dev = 20;
    filter_length = 50;
end
filter_length = std_dev*2;
x = [-filter_length:filter_length];
h = exp(-(x.^2)/(2*std_dev.^2))'; 
h = h';
%h = h*h';
h = h./sum(h(:));

symm_val = filter_length+1;
signal = symmetric_extensions(signal, [symm_val symm_val]);


%ping_signal = signal(l,:);
%filtered_ping_signal = filtfilt(h, 1, ping_signal);
%figure, plot(ping_signal); hold on
%plot(filtered_ping_signal, 'r');

filtered_signal = filter2(h, signal, 'same');
high_pass_signal = signal - filtered_signal;

%{
for l = 1:number_of_signals
    %l
    ping_signal = signal(l,:);
    filtered_ping_signal = filtfilt(h, 1, ping_signal);   
    high_pass_signal(l,:) = ping_signal - filtered_ping_signal;   
end
%}

filtered_signal = filtered_signal(symm_val + 1:end - symm_val,symm_val + 1:end - symm_val);
high_pass_signal = high_pass_signal(symm_val + 1:end - symm_val,symm_val + 1:end - symm_val);

end