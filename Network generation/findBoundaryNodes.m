%% FindPeriodicBond
% This code replaces a bind crossing the boundary of the box for two
% corresponding bonds in a systmatic and consistent way
% -----------------------------------------------------------------------
% Inputs:
%     a,b--> 3D coordinates of the periodic nodes;
%     dim--> Dimension of the problem

% Outputs
%     NewCoord--> Coordinates of the two extra nodes created
% =========================================================================
function NewCoords = FindPeriodicBond(a,b,dim)
    %% Assemble pulls
    NewCoords = zeros(2,3);
    pull = zeros(3,1);
    if (dim == 2)
        D = zeros(8,1);
        Directions = zeros(3,8);
    else 
        D = zeros(26,1);
        Directions = zeros(3,26); 
        lb = zeros(3,1);
        ub = zeros(3,1);        
    end

    if(dim == 2)
        M_pull = zeros(8,3);
        % Horizontal Displacements
        M_pull(1,1) = -1; M_pull(2,1) = 1;
        
        % Vertical displacemnts
        M_pull(3,2) = -1; M_pull(4,2) = 1;
        
        % Diagonal Displacements
        M_pull(5,1) = 1; M_pull(5,2) = 1; %(1,1)
        M_pull(6,1) = -1; M_pull(6,2) = 1; %(-1,1)
        M_pull(7,1) = -1; M_pull(7,2) = -1; %(-1,-1)
        M_pull(8,1) = 1; M_pull(8,2) = -1; %(1,-1)
    else
        M_pull = zeros(26,3);
        % Horizontal Displacements along x-axis
        M_pull(1,1) = -1; M_pull(2,1) = 1;
        
        % Vertical displacemnts along y-xis
        M_pull(3,2) = -1; M_pull(4,2) = 1;
        
        % Diagonal Displacements xy reference plane 
        M_pull(5,1) = 1; M_pull(5,2) = 1; %(1,1,0)
        M_pull(6,1) = -1; M_pull(6,2) = 1; %(-1,1,0)
        M_pull(7,1) = -1; M_pull(7,2) = -1; %(-1,-1,0)
        M_pull(8,1) = 1; M_pull(8,2) = -1; %(1,-1,0)

        % Diplacments along the z-axis
        M_pull(9,3) = -1; M_pull(10,3) = 1;

        % Horizontal displacemetns at "basement" (z = -1)
        M_pull(11,1) = -1; M_pull(11,3) = -1.; %(-1,0,-1)
        M_pull(12,1) = 1; M_pull(12,3) = -1.; %(1,0,-1)

        % Vertical Displacements at basement (z = -1)
        M_pull(13,2) = -1; M_pull(13,3) = -1.; %(0,-1,-1)
        M_pull(14,2) = 1;  M_pull(14,3) = -1.; %(0,1,-1)

        % Diagonal displacements on xy plane displaced to "basement"
        M_pull(15,1) = 1; M_pull(15,2) = 1; M_pull(15,3) = -1.; %(1,1,-1)
        M_pull(16,1) = -1; M_pull(16,2) = 1; M_pull(16,3) = -1.; %(-1,1,-1)
        M_pull(17,1) = -1; M_pull(17,2) = -1; M_pull(17,3) = -1.; %(-1,-1,-1)
        M_pull(18,1) = 1; M_pull(18,2) = -1; M_pull(18,3) = -1.; %(1,-1,-1)

        % Horizontal displacemetns at "ceilling" (z = 1)
        M_pull(19,1) = -1; M_pull(19,3) = 1.; %(-1,0,1)
        M_pull(20,1) = 1; M_pull(20,3) = 1.; %(1,0,1)

        % Vertical Displacements at "ceilling" (z = 1)
        M_pull(21,2) = -1; M_pull(21,3) = 1.; %(0,-1,1)
        M_pull(22,2) = 1;  M_pull(22,3) = 1.; %(0,1,1)

        % Diagonal displacements on xy plane displaced to "basement"
        M_pull(23,1) = 1; M_pull(23,2) = 1; M_pull(23,3) = 1.; %(1,1,1)
        M_pull(24,1) = -1; M_pull(24,2) = 1; M_pull(24,3) = 1.; %(-1,1,1)
        M_pull(25,1) = -1; M_pull(25,2) = -1; M_pull(25,3) = 1.; %(-1,-1,1)
        M_pull(26,1) = 1; M_pull(26,2) = -1; M_pull(26,3) = 1.; %(1,-1,1)

        % xyz limits of the original box
        ub = ones(3,1);

    end
    %% Mirror Nodes

    for i = 1:length(M_pull)
        M_b = b; pull = M_pull(i,:)';
        M_b = b + pull;

        v = M_b - a;
        D(i) = norm(v);
        Directions(:,i) = v;
    end
    %======================================================================

    %% Make Selection
    [d id] = min(D);

    v = Directions(:,id);
    
    pull = M_pull(id,:)';

    if(d>0.6)
        keyboard
    end

    % Set xyz axis of the selected box (3D only)
    if(dim == 3)
        lb = lb + pull; ub = ub + pull;
    end

    if(id == 1)
        % Box translated to (-1,0) is selected
        % Cross x = 0
        t = (- a(1))/(v(1));
        dummy1 = a + (t.*v);
        % Corresponding bonds crosses at x = 1
        dummy2 = dummy1;

    elseif(id == 2)
        % Box translated to (1,0) is selected
        % Cross x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        % Corresponding bonds crosses at x = 0
        dummy2 = dummy1;

    elseif(id == 3)
        % Box translated to (0,-1) is selected
        % Cross y = 0
        t = (- a(2))/(v(2));
        dummy1 = a + (t.*v);
        % Corresponding bonds crosses at y = 1
        dummy2 = dummy1;

    elseif(id == 4)
        % Box translated to (0,1) is selected
        % Cross y = 1
        t = (1. - a(2))/(v(2));
        dummy1 = a + (t.*v);
        % Corresponding bonds crosses at y = 0
        dummy2 = dummy1;

    elseif (id == 5)
        % Box translated by (1,1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = abs(dummy1) >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at y = 1
            t = (1. - a(2))/(v(2));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (1. - a(2))/(v(2));
            dummy2 = a + (t.*v);
            
        end
        
    elseif(id == 6)
        % Box translated by (-1,1) is selected
        % Assume ift crosses x = 0
        t = (- a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at y = 1
            t = (1. - a(2))/(v(2));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (1. - a(2))/(v(2));
            dummy2 = a + (t.*v);
            
        end
    elseif (id == 7)
        % Box tranaled by (-1,-1) is selected
        % Assume it crosses x = 0
        t = (-a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;
        
        if (length(find(check)) < length(dummy1))
            % it crosses at y = 0
            t = (-a(2))/(v(2));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (-a(2))/(v(2));
            dummy2 = a + (t.*v);
            
        end
    
    elseif (id == 8)
        % Box translated by (1,-1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at y = 0
            t = (- a(2))/(v(2));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (-a(2))/(v(2));
            dummy2 = a + (t.*v);
            
        end
    % ---------------------------------------------------------------------
    % This point foward is only for 3D
    elseif(id == 9)
        % Box traslated  by (0,0,-1) is selected        
        % Cross z = 0
        t = (- a(3))/(v(3));
        dummy1 = a + (t.*v);
        % Corresponding bonds crosses at z = 1
        dummy2 = dummy1;

    elseif(id == 10)
        % Box translated to (0,0,-1) is selected
        % Cross z = 1
        t = (1. - a(3))/(v(3));
        dummy1 = a + (t.*v);
        % Corresponding bonds crosses at z = 0
        dummy2 = dummy1;

    elseif(id == 11)
        % Box translated by (-1,0,-1) is selected
        % Assume ift crosses x = 0
        t = (-a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end
        
    elseif(id == 12)
        % Box translated by (1,0,-1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end

    elseif(id == 13)
        % Box translated by (0,-1,-1) is selected
        % Assume ift crosses y = 0
        t = (-a(2))/(v(2));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end

    elseif(id == 14)
        % Box translated by (0,1,-1) is selected
        % Assume ift crosses y = 1
        t = (1. - a(2))/(v(2));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end

    elseif(id == 15)
        % Box translated by (1,1,-1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 1 or z = 0.

                % Lets start with x = 1
                t = (1. - a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 0
                    t = (-a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 0, we need to now where the second point
            % is on the displaced box, either x = 1 or y = 1

            % Lets start with x = 1
            t = (1. - a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 0
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end

    elseif(id == 16)
        % Box translated by (-1,1,-1) is selected
        % Assume ift crosses x = 0
        t = (-a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 0 or z = 0.

                % Lets start with x = 0
                t = (-a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 0
                    t = (-a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 0, we need to now where the second point
            % is on the displaced box, either x = 0 or y = 1

            % Lets start with x = 0
            t = (-a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 0
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end

    elseif(id == 17)
        % Box translated by (-1,-1,-1) is selected
        % Assume ift crosses x = 0
        t = (-a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 0 or z = 0.

                % Lets start with x = 0
                t = (-a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 0
                    t = (-a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 0, we need to now where the second point
            % is on the displaced box, either x = 0 or y = 0

            % Lets start with x = 0
            t = (-a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 0
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end

    elseif(id == 18)
        % Box translated by (1,-1,-1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 0
            t = (-a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 0 or z = 0.

                % Lets start with x = 1
                t = (1. - a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 0
                    t = (-a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 0, we need to now where the second point
            % is on the displaced box, either x = 1 or y = 0

            % Lets start with x = 1
            t = (1. - a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 0
            t = (-a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end

    elseif(id == 19)
        % Box translated by (-1,0,1) is selected
        % Assume ift crosses x = 0
        t = (-a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end

    elseif(id == 20)
        % Box translated by (1,0,1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end        

     elseif(id == 21)
        % Box translated by (0,-1,1) is selected
        % Assume ift crosses y = 0 
        t = (-a(2))/(v(2));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end


    elseif(id == 22)
        % Box translated by (0,1,1) is selected
        % Assume ift crosses y = 1
        t = (1. - a(2))/(v(2));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = dummy1;
            dummy1 = a + (t.*v);
        else
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);    
        end

elseif(id == 23)
        % Box translated by (1,1,1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 1 or z = 1.

                % Lets start with x = 1
                t = (1. - a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 1
                    t = (1. - a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 1, we need to now where the second point
            % is on the displaced box, either x = 1 or y = 1

            % Lets start with x = 1
            t = (1. - a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end        

    elseif(id == 24)
        % Box translated by (-1,1,1) is selected
        % Assume ift crosses x = 0
        t = (-a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 0 or z = 1.

                % Lets start with x = 0
                t = (-a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 1
                    t = (1. - a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 1, we need to now where the second point
            % is on the displaced box, either x = 0 or y = 1

            % Lets start with x = 0
            t = (-a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 1
                t = (1. - a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end        

    elseif(id == 25)
        % Box translated by (-1,-1,1) is selected
        % Assume ift crosses x = 0
        t = (-a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 0 or z = 1.

                % Lets start with x = 0
                t = (-a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 1
                    t = (1. - a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 1, we need to now where the second point
            % is on the displaced box, either x = 0 or y = 0

            % Lets start with x = 0
            t = (-a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end

    elseif(id == 26)
        % Box translated by (1,-1,1) is selected
        % Assume ift crosses x = 1
        t = (1. - a(1))/(v(1));
        dummy1 = a + (t.*v);
        check = dummy1 >= 0. & dummy1 <= 1.0;

        if (length(find(check)) < length(dummy1))
            % check if it crosses at z = 1
            t = (1. - a(3))/(v(3));
            dummy1 = a + (t.*v);
            check = dummy1 >= 0. & dummy1 <= 1.0;

            if(length(find(check)) < length(dummy1))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy1 = a + (t.*v);
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % If it crosses at y = 1, we need to now where the second point
                % is on the displaced box, either x = 0 or z = 1.

                % Lets start with x = 1
                t = (1. - a(1))/(v(1));
                dummy2 = a + (t.*v);

                check = dummy2 >= lb & dummy2 <= ub;

                if(length(find(check)) < length(dummy2))
                    % It crosses at z = 1
                    t = (1. - a(3))/(v(3));
                    dummy2 = a + (t.*v);
                end
            end

            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % If it crosses at z = 1, we need to now where the second point
            % is on the displaced box, either x = 1 or y = 0

            % Lets start with x = 1
            t = (1. - a(1))/(v(1));
            dummy2 = a + (t.*v);

            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        else
            % Our guess was right and we need to check which other plane
            % will be intercepted on the correspndinf box

            % Try z = 1
            t = (1. - a(3))/(v(3));
            dummy2 = a + (t.*v);
            
            check = dummy2 >= lb & dummy2 <= ub;

            if(length(find(check)) < length(dummy2))
                % It crosses at y = 0
                t = (-a(2))/(v(2));
                dummy2 = a + (t.*v);
            end

        end

    end

    

    %% Compute the corresponding points
    % Translate dummy2 to the original box 
    d_loss = norm(dummy2 - dummy1); % Usually is zero
    dummy2 = dummy2 - pull;

    % Check for rounding errors
    for k = 1:dim 
        % Check for the first coordinate
        if (dummy1(k) < 0. || dummy1(k) > 1.)
            if(abs(dummy1(k)) < 1e-6)
                dummy1(k) = 0.;
            else
                dummy1(k) = 1.;
            end
        end
        
        % Second point
        if (dummy2(k) < 0. || dummy2(k) > 1.)
            if(abs(dummy2(k)) < 1e-6)
                dummy2(k) = 0.;
            else
                dummy2(k) = 1.;
            end
        end
    end

    NewCoords(1,:) = dummy1;
    NewCoords(2,:) = dummy2;

    


end