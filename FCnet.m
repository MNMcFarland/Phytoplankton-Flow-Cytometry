classdef FCnet
%
% A simple feedforward artificial neural network (multilayer perceptron)
% for classifying flow cytometry data (events x parameters, i.e. inputs as row vectors)
% with optional dropout specified for each internal layer
%
% USAGE: 
%   net = FCnet(sizes); % create a new network with random initialization
%   net = FCnet(sizes, drprb); % create a new network with dropout
%   net = train(net, train_in, train_lb); % train the network 
%   lbl = predict(net, X); % classify data in X with trained network
%
% INPUT
%   sizes - size of each layer (vector), sizes(1) = size of input (# of parameters)
%   drprb - (optional) drop probabilities for each hidden layer (vector, same size as sizes)
%           all must satisfy 0 <= drprb < 1, ignored for 1st and last layer
%
% OUTPUT
%   net - untrained network structure with built in functions:
%       predict(net, input) - use the network to classify data in 'input'
%       train(net, train_in, train_lb) - train the network
%
    properties
        nL % = length(sizes); % number of layers
        sizes 
        biases
        weights
        drp % = zeros(size(sizes)); % drop probabilities for each layer
        % hyperparameter defaults
        epochs = 50; % default number of epochs (passes through data set during training)
        mbs = 50; % default mini batch size
        eta = .001; % default learning rate
        lmbda = 0; % length(train_labels)*1e-4; % regularization parameter, 0 = no regularization (default)
        n % size of training data set
        lblID
        lblnames = {};
        accuracy
        % mu = 0.5; % velocity friction parameter (momentum coefficient) range = 0:1, 0 = no velocity effect
        %  epoch % for Adam optimizer
    end

    methods (Static)
        function a = act(z)
%             a = max(0,z); % Rectified linear unit (ReLU) activation function
            % a = 1./(1+exp(-z)); % logistic (sigmoid) activation function
            a = z./(1+exp(-z)); % swish activation function
        end
        %%
        function a = act_prime(z)
            % derivative of the ReLU activation function
%             a = zeros(size(z));
%             a(z>0) = 1;
            
            % a = z.*(1-z); % derivative of the sigmoid activation function

            % derivative of the swish activation function
            z = z./(1+exp(-z));
            a = z + (1-z)./(1+exp(-z)); % derivative of swish activation function
        end
        %%
        function s = sftmx(z)
            % numerically stabilized softmax
            z = z-max(z,[],2);
            ez = exp(z);
            s = ez./sum(ez,2,'omitnan');
        end
    end
    methods
        function obj = FCnet(sz,drp) % constructor
            if nargin > 0
                obj.sizes = sz;
                obj.nL = length(sz);
                if ~exist('drp','var') || isempty(drp) 
                    obj.drp = zeros(size(obj.sizes));
                elseif length(drp)==length(obj.sizes)
                    obj.drp = drp;
                else
                    error('FCnet_class input error')
                end
                obj = init(obj);
            end
        end
        %%
        function obj = init(obj)
            if isempty(obj.nL); obj.nL = length(obj.sizes); end
            % initial random weights and biases
            for l = 2:length(obj.sizes)
            %     obj.biases{l} = randn(1,obj.sizes(l)); % biases for each neuron in each layer
                obj.biases{l} = zeros(1,obj.sizes(l)); % biases for each neuron in each layer
                obj.weights{l} = randn(obj.sizes(l-1), obj.sizes(l)) * sqrt(2/obj.sizes(l-1)); % weights for each input to each neuron (neuron x input) in each layer 
            %     obj.weights{l} = randn(obj.sizes(l-1), obj.sizes(l))/sqrt(obj.sizes(l-1)); % weights for each input to each neuron (neuron x input) in each layer 
            %     obj.biases{l} = rand(1,sizes(l))'-.5; % biases for each neuron in each layer
            %     obj.weights{l} = rand(sizes(l), sizes(l-1))-.5; % weights for each input to each neuron (neuron x input) in each layer 
            end
        end
        %%
        function [lbl,prb] = predict(obj,a)
            % Classify input by applying weights, biases, and activation function to each layer
            % of the network (excluding input layer)
            %
            % INPUT
            %   obj - network object
            %   a - input with observations (events) as row vectors (events x parameters)
            %
            % OUTPUT
            %   lbl - predicted class labels
            %   a - softmax activations for each input and class
            %   prb - probability of predictions
            
            for l = 2:obj.nL-1
                a = obj.act(a*obj.weights{l} + obj.biases{l});
            end
            a = obj.sftmx(a*obj.weights{end} + obj.biases{end}); % softmax final layer
            
            [prb,lbl] = max(a,[],2); % class labels and probabilities as vectors
            lbl = obj.lblID(lbl);
        end
        %%
        function obj = train(obj, train_in, train_lb, test_in, test_lb, epochs, mbs, eta, lmbda)
            % Train network using Stochastic Gradient Descent
            %
            % INPUT:
            %   obj - network object
            %   train_in - input data as row vectors
            %   train_lb - label vector with ID for each row of train_in
            %   test_in - data to evaluate accuracy
            %   test_lb - labels to evaluate accuracy
            %   epochs - the number of epochs to train for
            %   mbs - size of the mini-batches to use when sub-sampling
            %   eta - learning rate
            %   lmbda - L2 regularization parameter
            %
            % OUTPUT:
            %   obj - trained network object
            %   result - accuracy determined from test input
            
            obj.n = length(train_lb); % size of training data set 
            
            if exist('epochs','var') && ~isempty(epochs); obj.epochs = epochs; end
            if exist('mbs','var') && ~isempty(mbs); obj.mbs = mbs; end
            if exist('eta','var') && ~isempty(eta); obj.eta = eta; end
            if exist('lmbda','var') && ~isempty(lmbda); obj.lmbda = lmbda; end
            
            if ~exist('test_in','var'); test_in = []; end
            if ~exist('test_lb','var'); test_lb = []; end
            
            % desired network output (convert to one hot encoding)
            train_out = zeros(obj.n,obj.sizes(end));
            % idx = sub2ind(size(train_out),(1:n)',train_lb);
            if isempty(obj.lblID) % generate new label IDs from training data
                [obj.lblID,~,rank_lb] = unique(train_lb);
            else % use existing label IDs (for re-training, make sure IDs match existing classes!)
                [~,~,rank_lb] = unique(train_lb);
            end
            idx = sub2ind(size(train_out),(1:obj.n)',rank_lb);
            train_out(idx) = 1;
            
            % determine weights to correct for imbalance
            % clswts = sum(train_out)/obj.n; % class weights
            % evtwts = clswts(train_lb)'; % class weight for each event
            % evtwts = 1 - evtwts.*train_out; % weight array for each output neuron
            
            % initialize velocity parameters for momentum based gradient descent
            % for l=2:obj.nL
            %     obj.vw{l}=zeros(size(obj.weights{l},1),1);
            %     obj.vb{l}=zeros(size(obj.biases{l},1),1);
            % end
            
            % initialize Adam parameters - doesn't work yet
            % obj.beta1 = 0.9;
            % obj.beta2 = 0.999;
            % obj.epsln = 10^-8;
            % obj.m = num2cell(zeros(1,obj.nL)); % cell(1,obj.nL);
            % obj.v = num2cell(zeros(1,obj.nL)); % cell(1,obj.nL);
            
            obj.accuracy = zeros(1,obj.epochs);
            for j = 1:obj.epochs
                %   obj.epoch = j;
                idx = randperm(obj.n)'; % shuffled index
                mb_in = train_in(idx,:); % shuffle train data
                mb_lb = train_out(idx,:); % shuffle train labels
                for k = obj.mbs:obj.mbs:obj.n
                    obj = backprop(obj, mb_in(k-mbs+1:k,:), mb_lb(k-mbs+1:k,:));
                %     obj = backprop(obj, mb_in(k-obj.mbs+1:k,:), mb_lb(k-obj.mbs+1:k,:), evtwts(k-obj.mbs+1:k,:)); 
                end
                disp(['Epoch ' num2str(j) ' training complete'])
                if ~isempty(test_in) && ~isempty(test_lb)
                    obj.accuracy(j) = evaluate(obj,test_in,test_lb);
                    disp(['accuracy = ' num2str(obj.accuracy(j)) '%'])
%                     obj.accuracy(j) = result(j);
                end
            end
        end
        %%
        function obj = backprop(obj, in, lb, evtwts)
            % Backpropagation vectorized over mini batch of training data
            %
            % INPUT:
            %   obj - network object containing weights and biases to be updated
            %   in - mini batch of network input, events x parameters 
            %   lb - one-hot encoded known labels, desired network output 
            %   evtwts - class weights for each event
            %
            % OUTPUT:
            %   obj - network object with updated weights and biases
            
            if ~exist('evtwts','var') || isempty(evtwts); evtwts = 1; end
            
            % feed forward while remembering activations and zs for each layer
            zs = cell(1,obj.nL);
            activations = cell(1,obj.nL);
            activations{1} = in;
            msk = cell(1,obj.nL);
            for l = 2:obj.nL-1
                zs{l} = activations{l-1} * obj.weights{l} + obj.biases{l};
                activations{l} = obj.act(zs{l});
                if obj.drp(l) > 0 && obj.drp(l) < 1
                    msk{l} = binornd(1, 1-obj.drp(l), size(activations{l})) / (1-obj.drp(l)); % dropout mask
                    activations{l} = activations{l} .* msk{l}; % apply dropout
                end
            end
            zs{end} = activations{l} * obj.weights{end} + obj.biases{end};
            activations{end} = obj.sftmx(zs{end}); % softmax final layer
            
            % Loss and Cost
            % L = vecnorm(lb-activations{end},2,2)^2/2; % quadratic loss
            % [~,id] = max(lb,[],2); % column index of one-hot labels = class number
            % lid = sub2ind(size(activations{end}),1:obj.mbs,id);
            % L = -log(activations{end}(lid)+eps); % cross entropy loss = -log likelihood loss
            % % L = -sum(lb.*log(activations{end}+eps),2); % cross entropy loss = -log likelihood loss
            % C = sum(L)/obj.mbs; % Cost
            
            % backward propagation to determine mean gradients
            % mbs = size(lb,1); % mini batch size to compute mean over mini batch
            % dC = activations{end}-lb; % derivative of the cost function
            % delta = activations{end}-lb; % delta for cross entropy cost with sigmoid activation function = dC/dz{L}
            % delta = (activations{end}-lb).*obj.sftmx_prime(zs{end}); % delta for -log likelihood cost with softmax activation function = dC/dz{L}
            delta = activations{end} - lb; % delta for -log likelihood cost with softmax activation function = dC/dz{L}
            
            delta = delta.*evtwts; % apply weights to account for imbalanced training data
            
            nabla_b{obj.nL} = sum(delta,1);%/mbs;
            nabla_w{obj.nL} = activations{obj.nL-1}'*delta;%/mbs;
            for l = obj.nL-1:-1:2
                delta = (delta*obj.weights{l+1}').*obj.act_prime(zs{l});
                if obj.drp(l) > 0 && obj.drp(l) < 1
                    delta = delta .* msk{l}; % dropout
                end
                nabla_b{l} = sum(delta,1);%/mbs;
                nabla_w{l} = activations{l-1}'*delta;%/mbs;
            end
            
            % update weights and biases
            for l = 2:obj.nL
            %     obj.weights{l} = obj.weights{l} - (eta/mbs) * nabla_w{l}; % no regularization
                obj.weights{l} = (1-obj.eta*(obj.lmbda/obj.n)) * obj.weights{l} - obj.eta / obj.mbs * nabla_w{l}; % with L2 regularization
                obj.biases{l} = obj.biases{l} - obj.eta / obj.mbs * nabla_b{l};
                
                % momentum based gradient descent
            %     obj.vw{l} = obj.mu*obj.vw{l} - (eta/m) * nabla_w{l}; % velocity
            %     obj.vb{l} = obj.mu*obj.vb{l} - (eta/m) * nabla_b{l}; % velocity
            %     obj.weights{l} = (1-eta*(lmbda/n))*obj.weights{l} + obj.vw{l}; % with L2 regularization and velocity
            %     obj.biases{l} = obj.biases{l} + obj.vb{l};
            
                % Adam - doesn't work yet
            %     obj.m{l} = obj.beta1*obj.m{l} + (1-obj.beta1)*nabla_w{l};
            %     obj.v{l} = obj.beta2*obj.v{l} + (1-obj.beta2)*nabla_w{l}.^2;
            %     obj.m{l} = obj.m{l}/(1-obj.beta1^obj.epoch);
            %     obj.v{l} = obj.v{l}/(1-obj.beta2^obj.epoch);
            % %     obj.weights{l} = obj.weights{l} - obj.eta * obj.m{l} ./ (sqrt(obj.v{l})+obj.epsln);
            %     obj.weights{l} = obj.weights{l} - obj.eta / obj.mbs * obj.m{l} ./ (sqrt(obj.v{l})+obj.epsln);
            end
        end
        %%
        function result = evaluate(obj, test_in, test_lb)
            % evaluate a batch of input and determine accuracy
            % [~,test_results] = max(net.predict(net,test_in),[],2);
            [test_results] = predict(obj, test_in);
            result = sum(test_results==test_lb)/length(test_lb)*100;
        end
    end
end