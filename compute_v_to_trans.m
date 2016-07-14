function [ r_ax,r_ang,scale ] = compute_v_to_trans( Vs )
scale=sqrt(sum(Vs.^2,2));
Vs=[Vs(:,1)./scale Vs(:,2)./scale Vs(:,3)./scale];
scale=scale/2;
ax=[0 1 0];
axs=repmat(ax,size(Vs,1),1);
 crs=cross(Vs,axs);
 %first arg should be size of cross-prod but all Vs were normalized
 r_ang=-atan2(sqrt(sum(crs.^2,2)),dot(Vs',axs')');
 r_ax=normr(crs);
end

