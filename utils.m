% utils.m
function varargout = utils(func, varargin)
    switch func
        case 'buildAgent'
            varargout{1} = buildAgent(varargin{:});
        case 'iofl'
            varargout{1} = iofl(varargin{:});
        otherwise
            error('utils: unknown function "%s"', func);
    end
end

% -------------------------------------------------------------------------

function agent = buildAgent(phi, B)
    syms p_x p_y v_x v_y u_x u_y real
    p = [p_x; p_y];
    v = [v_x; v_y];
    u = [u_x; u_y];
    x = [p_x; v_x; p_y; v_y];

    phi_var = phi(p, v);
    B_var   = B(p, v);
    xdot    = phi_var + B_var * u;
    y       = [p_x; p_y];

    agent.phi  = phi_var;
    agent.B    = B_var;
    agent.x    = x;
    agent.xdot = matlabFunction(xdot, 'Vars', {x, u});
    agent.y    = simplify(y);
end

% -------------------------------------------------------------------------

function u = iofl(agent)
    agent_name  = inputname(1);
    phi         = agent.phi;
    B           = agent.B;
    y           = agent.y;
    x           = agent.x;

    num_outputs = length(y);
    num_inputs  = size(B, 2);
    max_rel_deg = length(phi);

    beta  = sym(zeros(num_outputs, num_inputs));
    alpha = sym(zeros(num_outputs, 1));
    r     = zeros(num_outputs, 1);

    for idx_out = 1:num_outputs
        h_curr        = y(idx_out);
        rel_deg_found = false;

        for ga = 1:max_rel_deg
            Lg_h = jacobian(h_curr, x) * B;

            if any(~isAlways(Lg_h == 0, 'Unknown', 'false'))
                beta(idx_out, :) = Lg_h;
                alpha(idx_out)   = jacobian(h_curr, x) * phi;
                r(idx_out)       = ga;
                rel_deg_found    = true;
                break;
            end

            h_curr = jacobian(h_curr, x) * phi;
        end

        if ~rel_deg_found
            error('FeedbackLinearization:UndefinedRelativeDegree', ...
                'Relative degree undefined for output %d within %d iterations.', ...
                idx_out, max_rel_deg);
        end
    end

    rnk = rank(beta);
    if rnk < num_outputs
        error('FeedbackLinearization:SingularDecouplingMatrix', ...
            'Decoupling matrix is singular.');
    end

    syms nu_x nu_y
    nu    = [nu_x; nu_y];
    u_sym = simplify(beta \ (nu - alpha));
    u     = matlabFunction(u_sym, 'Vars', {x, nu});

    fprintf('\n/--- IOFL Summary for agent %s ---\n\n', agent_name);
    fprintf('Vector relative degree: r = [%s]\n\n', num2str(r.'));
    fprintf('Drift term (alpha):\n');
    alpha_str = evalc('disp(vpa(alpha, 3))');
    printIndented(alpha_str, '    ');
    fprintf('Decoupling matrix (beta):\n');
    beta_str = evalc('disp(vpa(beta, 3))');
    printIndented(beta_str, '    ');
    fprintf('Rank: %d ( == dim(u) ). Agent is linearizable.\n\n', rnk);
    fprintf('Control law transformation:\nu = \n');
    u_string = string(vpa(u_sym, 3));
    for idx_out = 1:num_outputs
        fprintf('    %s\n', u_string(idx_out));
    end
    fprintf('\n\\---------------------------------\n\n\n');
end

% -------------------------------------------------------------------------

function printIndented(raw_text, indentation)
    lines = splitlines(strtrim(raw_text));
    for i = 1:length(lines)
        if ~isempty(lines{i})
            fprintf('%s%s\n', indentation, lines{i});
        end
    end
    fprintf('\n');
end