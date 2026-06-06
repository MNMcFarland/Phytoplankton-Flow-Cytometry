function net = FCnet_gpu(sizes,drprb)
%%
% A simple feedforward artificial neural network (multilayer perceptron)
% for classifying flow cytometry data (events x parameters, i.e. inputs as row vectors)
% with dropout
%
% USAGE: 
%   net = FCnet(sizes); % create a new network
%   [net, result] = net.SGD(net, train_in, train_lb); % train the network 
%   lbl = net.feedforward(net, X); % classify data in X with trained network
%
% INPUT
%   sizes - size of each layer (vector), sizes(1) = size of input (# of parameters)
%   drprb - (optional) drop probabilities for each hidden layer (vector, same size as sizes)
%           all must satisfy 0 <= drprb < 1, ignored for 1st and last layer
%
% OUTPUT
%   net - untrained network structure with built in functions:
%       net.feedforward(net, input) - run the network with 'input'
%       net.SGD(net, train_in, train_lb) - train the network
%

net.nL = gpuArray(length(sizes)); % number of layers
net.sizes = gpuArray(sizes);

if exist('drprb','var') && ~isempty(drprb)
    net.drp = gpuArray(drprb); % dropout probabilities
else
    net.drp = zeros(size(sizes),'gpuArray');
end

% initial random weights and biases
for l = 2:length(net.sizes)
%     net.biases{l} = randn(1,net.sizes(l)); % biases for each neuron in each layer
    net.biases{l} = zeros(1,net.sizes(l)); % biases for each neuron in each layer
    net.weights{l} = randn(net.sizes(l-1), net.sizes(l)) * sqrt(2/net.sizes(l-1)); % weights for each input to each neuron (neuron x input) in each layer 
%     net.weights{l} = randn(net.sizes(l-1), net.sizes(l))/sqrt(net.sizes(l-1)); % weights for each input to each neuron (neuron x input) in each layer 
%     net.biases{l} = rand(1,sizes(l))'-.5; % biases for each neuron in each layer
%     net.weights{l} = rand(sizes(l), sizes(l-1))-.5; % weights for each input to each neuron (neuron x input) in each layer 
    net.biases{l} = gpuArray(net.biases{l});
    net.weights{l} = gpuArray(net.weights{l});
end

% hyperparameter defaults
net.epochs = gpuArray(50); % default number of epochs (passes through data set during training)
net.mbs = gpuArray(50); % default mini batch size
net.eta = gpuArray(.001); % default learning rate
net.lmbda = gpuArray(0); % length(train_labels)*1e-4; % regularization parameter, 0 = no regularization (default)
% net.mu = 0.5; % velocity friction parameter (momentum coefficient) range = 0:1, 0 = no velocity effect

% activation function
% net.act = @(z) 1./(1+exp(-z)); % logistic (sigmoid) activation function
% net.act_prime = @(z) net.act(z).*(1-net.act(z)); % derivative of sigmoid activation function
net.act = @(z) max(0,z); % Rectified linear unit (ReLU) activation function
net.act_prime = @ReLU_prime; % derivative of ReLU activation function
% net.act = @swish; % swish activation function
% net.act_prime = @swish_prime; % derivative of swish activation function
net.sftmx = @sftmx;
% net.sftmx_prime = @sftmx_prime;

% function handles
net.feedforward = @feedforward;
net.SGD = @SGD;
net.update_mini_batch = @update_mini_batch;
net.evaluate = @evaluate;

function a = ReLU_prime(z)
% derivative of the ReLU activation function
a = zeros(size(z),'gpuArray');
a(z>0) = 1;

function s = sftmx(z)
% numerically stabilized softmax
z = z-max(z,[],2);
ez = exp(z);
s = ez./nansum(ez,2);

% function ds = sftmx_prime(z) % not actually needed if only used at output layer
% z = z-max(z,[],2);
% ez = exp(z);
% s = ez./nansum(ez,2);
% ds1 = s(1,:)*(eye(length(s(1,:)))-s(1,:)'); % this doesn't work for matrix input, row vector only
% sr = s(:); % reshape
% ds = diag(sr) - sr*sr'; % returns the Jacobian, what to do with this?

% function a = swish(z)
% a = z./(1+exp(-z)); % swish activation function
% 
% function a = swish_prime(z)
% z = z./(1+exp(-z));
% a = z + (1-z)./(1+exp(-z)); % derivative of swish activation function

function [a,lbl,prb] = feedforward(net,a)
%% Classify input by applying acitivation function to each layer of the network (excluding input layer)
%
% INPUT
%   net - network structure
%   a - input with observations (events) as row vectors (events x parameters)
%
% OUTPUT
%   a - softmax activations for each input
%   lbl - predicted class labels
%   prb - probability of predictions

a = gpuArray(a);
for l = 2:net.nL-1
    a = net.act(a*net.weights{l} + net.biases{l});
end
a = net.sftmx(a*net.weights{end} + net.biases{end}); % softmax final layer

[prb,lbl] = max(a,[],2); % class labels and probabilities as vectors
lbl = net.lblnames(lbl);

function [net,result] = SGD(net, train_in, train_lb, test_in, test_lb, epochs, mbs, eta, lmbda)
%% Train network using Stochastic Gradient Descent
% INPUT
%   net - network structure
%   train_in - input data as row vectors
%   train_lb - label vector for each row of train_in
%   test_in - data to evaluate accuracy
%   test_lb - labels to evaluate accuracy
%   epochs - the number of epochs to train for
%   mbs - size of the mini-batches to use when sub-sampling
%   eta - learning rate
%   lmbda - L2 regularization parameter
%
% OUTPUT
%   net - trained network structure
%   result - accuracy determined from test input

train_in = gpuArray(train_in);
train_lb = gpuArray(train_lb);
test_in = gpuArray(test_in);
test_lb = gpuArray(test_lb);

n = length(train_lb);
net.n = gpuArray(n); % size of training data set 

if exist('epochs','var') && ~isempty(epochs); net.epochs = gpuArray(epochs); end
if exist('mbs','var') && ~isempty(mbs); net.mbs = gpuArray(mbs); end
if exist('eta','var') && ~isempty(eta); net.eta = gpuArray(eta); end
% if ~exist('lmbda','var') || isempty(lmbda); end
if exist('lmbda','var') && ~isempty(lmbda); net.lmbda = gpuArray(lmbda); end

if ~exist('test_in','var'); test_in = []; end
if ~exist('test_lb','var'); test_lb = []; end

% desired network output (convert to one hot encoding)
train_out = zeros(n,net.sizes(end),'gpuArray');
% idx = sub2ind(size(train_out),(1:n)',train_lb);
[net.lblnames,~,rank_lb] = unique(train_lb);
net.lblnames = gpuArray(net.lblnames);
idx = sub2ind(size(train_out),(1:n)',rank_lb);
train_out(idx) = 1;

% determine class weights to correct for imbalance
% clswts = sum(train_out)/net.n; % class weights
% evtwts = clswts(train_lb)'; % class weight for each event
% evtwts = 1 - evtwts.*train_out; % weight array for each output neuron

% initialize velocity parameters for momentum based gradient descent
% for l=2:net.nL
%     net.vw{l}=zeros(size(net.weights{l},1),1);
%     net.vb{l}=zeros(size(net.biases{l},1),1);
% end

% initialize Adam parameters - doesn't work yet
% net.beta1 = 0.9;
% net.beta2 = 0.999;
% net.epsln = 10^-8;
% net.m = num2cell(zeros(1,net.nL)); % cell(1,net.nL);
% net.v = num2cell(zeros(1,net.nL)); % cell(1,net.nL);

result = zeros(1,net.epochs,'gpuArray');
for j = 1:net.epochs
    net.epoch = gpuArray(j);
    idx = gpuArray(randperm(net.n)'); % shuffled index
    mb_in = train_in(idx,:); % shuffle train data
    mb_lb = train_out(idx,:); % shuffle train labels
    for k = net.mbs:net.mbs:net.n
        net = update_mini_batch(net, mb_in(k-mbs+1:k,:), mb_lb(k-mbs+1:k,:));
%         net = update_mini_batch(net, mb_in(k-net.mbs+1:k,:), mb_lb(k-net.mbs+1:k,:), evtwts(k-net.mbs+1:k,:)); 
    end
    disp(['Epoch ' num2str(j) ' training complete'])
    if ~isempty(test_in) && ~isempty(test_lb)
        result(j) = evaluate(net,test_in,test_lb);
        disp(['accuracy = ' num2str(result(j)) '%'])
        net.accuracy = result(j);
    end
end

function net = update_mini_batch(net, in, lb, evtwts)
%% Backpropagation vectorized over mini batch of training data
% INPUT:
%   net - network structure containing weights and biases to be updated
%   in - mini batch of network input, events x parameters 
%   lb - one-hot encoded known labels, desired network output 
%   evtwts - class weights for each event
%
% OUTPUT:
%   net - network structure with updated weights and biases

if ~exist('evtwts','var') || isempty(evtwts); evtwts = gpuArray(1); end

% feed forward while remembering activations and zs for each layer
n = gather(net.nL);
zs = cell(1,n);
activations = cell(1,n);
activations{1} = in;
msk = cell(1,n);
for l = 2:n-1
    zs{l} = activations{l-1} * net.weights{l} + net.biases{l};
    activations{l} = net.act(zs{l});
    if net.drp(l) > 0 && net.drp(l) < 1
        msk{l} = binornd(1, 1-net.drp(l), size(activations{l})) / (1-net.drp(l)); % dropout mask
        msk{l} = gpuArray(msk{l});
        activations{l} = activations{l} .* msk{l}; % apply dropout
    end
end
zs{end} = activations{l} * net.weights{end} + net.biases{end};
activations{end} = net.sftmx(zs{end}); % softmax final layer

% Loss and Cost
% L = vecnorm(lb-activations{end},2,2)^2/2; % quadratic loss
% [~,id] = max(lb,[],2); % column index of one-hot labels = class number
% lid = sub2ind(size(activations{end}),1:net.mbs,id);
% L = -log(activations{end}(lid)+eps); % cross entropy loss = -log likelihood loss
% % L = -sum(lb.*log(activations{end}+eps),2); % cross entropy loss = -log likelihood loss
% C = sum(L)/net.mbs; % Cost

% backward propagation to determine mean gradients
% mbs = size(lb,1); % mini batch size to compute mean over mini batch
% dC = activations{end}-lb; % derivative of the cost function
% delta = activations{end}-lb; % delta for cross entropy cost with sigmoid activation function = dC/dz{L}
delta = activations{end} - lb; % delta for -log likelihood cost with softmax activation function = dC/dz{L}
% delta = (activations{end}-lb).*net.sftmx_prime(zs{end}); % delta for -log likelihood cost with softmax activation function = dC/dz{L}

delta = delta.*evtwts; % apply weights to account for imbalanced training data

nabla_b{net.nL} = sum(delta,1);%/mbs;
nabla_w{net.nL} = activations{net.nL-1}'*delta;%/mbs;
for l = net.nL-1:-1:2
    delta = (delta*net.weights{l+1}').*net.act_prime(zs{l});
    if net.drp(l) > 0 && net.drp(l) < 1
        delta = delta .* msk{l}; % dropout
    end
    nabla_b{l} = sum(delta,1);%/mbs;
    nabla_w{l} = activations{l-1}'*delta;%/mbs;
end

% update weights and biases
for l = 2:n
%     net.weights{l} = net.weights{l} - (eta/mbs) * nabla_w{l}; % no regularization
    net.weights{l} = (1-net.eta*(net.lmbda/net.n)) * net.weights{l} - net.eta / net.mbs * nabla_w{l}; % with L2 regularization
    net.biases{l} = net.biases{l} - net.eta / net.mbs * nabla_b{l};
    
    % momentum based gradient descent
%     net.vw{l} = net.mu*net.vw{l} - (eta/m) * nabla_w{l}; % velocity
%     net.vb{l} = net.mu*net.vb{l} - (eta/m) * nabla_b{l}; % velocity
%     net.weights{l} = (1-eta*(lmbda/n))*net.weights{l} + net.vw{l}; % with L2 regularization and velocity
%     net.biases{l} = net.biases{l} + net.vb{l};

    % Adam - doesn't work yet
%     net.m{l} = net.beta1*net.m{l} + (1-net.beta1)*nabla_w{l};
%     net.v{l} = net.beta2*net.v{l} + (1-net.beta2)*nabla_w{l}.^2;
%     net.m{l} = net.m{l}/(1-net.beta1^net.epoch);
%     net.v{l} = net.v{l}/(1-net.beta2^net.epoch);
% %     net.weights{l} = net.weights{l} - net.eta * net.m{l} ./ (sqrt(net.v{l})+net.epsln);
%     net.weights{l} = net.weights{l} - net.eta / net.mbs * net.m{l} ./ (sqrt(net.v{l})+net.epsln);
end

function result = evaluate(net, test_in, test_lb)
%% evaluate a batch of input and determine accuracy
% [~,test_results] = max(net.feedforward(net,test_in),[],2);
[~,test_results] = net.feedforward(net,test_in);
result = sum(test_results==test_lb)/length(test_lb)*100;


