% Created by Isabel Llorente-Garcia, 2011.

% Create circular mask around position of found and accepted spots:
        % circle with ones inside and zeros outside. Default radius is 5 pixels.
        % (circle_mask is a square image of the same size as frameROI with
        % a circular mask around the position of the spot):
        circle_mask = zeros(size(frameROI,1),size(frameROI,2)); % pre-define circle_mask size, initialise.       
        % Need error control when taking points at edge: do not use points
        % too close to edge:
        if (round(x_centre)-2*params.inner_radius)<1 || (round(x_centre)+2*params.inner_radius)>size(frameROI,2) || (round(y_centre)-2*params.inner_radius)<1 || (round(y_centre)+2*params.inner_radius)>size(frameROI,1)
            % error control for when spots are too close to edge of image ROI.
            % Do nothing (skip spot) if too close to edge.
        else
            % only for a small square (saves time) of size 2*params.inner_radius+1 around the
            % found spot, build circular mask:
            for ii = (round(y_centre)-2*params.inner_radius):(round(y_centre)+2*params.inner_radius) % ii = index for rows
                for jj = (round(x_centre)-2*params.inner_radius):(round(x_centre)+2*params.inner_radius) % jj = index for columns
                    % circle mask: if distance to point (x_centre,y_centre) is <=params.inner_radius, then 1, otherwise, 0:
                    sum_sq = (jj-x_centre)^2 + (ii-y_centre)^2;
                    if round(sqrt(sum_sq))<=params.inner_radius
                        circle_mask(ii,jj)=1;
                    else
                        circle_mask(ii,jj)=0;
                    end
                end
            end
            nnz(circle_mask) % number of non-zero elements in matrix circle_mask.
%                      disp(k)
%                      disp(n)
%                      size(frameROI)
%                      size(circle_mask)
            %         figure;
            %         imshow(circle_mask,[],'Border','tight')
            %         figure;
            %         imshow(circle_mask.*frameROI,[],'Border','tight')
            
            Ispots(k,n) = sum(sum(circle_mask.*frameROI));            
        end        