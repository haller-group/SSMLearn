function RBF_object = RBF_interpolator(X, y)
        % given a set of training data such that f(X) = y, we compute the
        % radial basis function interpolation of f with a linear kernel
       kernel_matrix = -pdist2(X', X'); % X is size (n_inputs, n_samples)
       coeffs = kernel_matrix \ y';
       RBF_object.coeffs = coeffs;
       RBF_object.centers = X;
       RBF_object.kernel_querry = @(x_query) -pdist2(x_query', RBF_object.centers');
       RBF_object.evaluate = @(x_query) transpose(RBF_object.kernel_querry(x_query) * RBF_object.coeffs);
end