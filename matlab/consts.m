%%
if ~exist('kyber', 'var')
    kyberk = 3;
    kyber = kyberk*256;
end

%% roundoff errors

switch kyber
    case 512
        err_round_v = dlmread('data/err_round_v_512.txt');
    case 768
        err_round_v = dlmread('data/err_round_v_768.txt');
    case 1024
        err_round_v = dlmread('data/err_round_v_1024.txt');
    otherwise
        fprintf('unknown parameter set!\n');
end
err_round_v(:, 2) = err_round_v(:, 2)/sum(err_round_v(:, 2));

% err_round_u is so small, we can hardcode it in this script
switch kyber
    case {512, 768}
        err_round_u = [
        -2 128
        -1 1024
        0 1024
        1 1024
        2 129
        ];
    case 1024
        err_round_u = [
        -1 640
        0 2048
        1 641
        ];
    otherwise
        fprintf('unknown parameter set!\n');
end

err_round_u(:, 2) = err_round_u(:, 2)/sum(err_round_u(:, 2));

%% empirical error

switch kyber
    case 512
        err_empirical = dlmread('data/err_empirical_512.txt');
    case 768
        err_empirical = dlmread('data/err_empirical_768.txt');
    case 1024
        err_empirical = dlmread('data/err_empirical_1024.txt');
    otherwise
        fprintf('unknown parameter set!\n');
end
err_empirical(:, 2) = err_empirical(:, 2)/sum(err_empirical(:, 2));
