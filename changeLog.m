%% Change log for version saved on 2016/12/1
% This version has two changes: the use of intensity in the particle
% potential well and the modification of the 1/r part of the interparticle
% potential interaction.
 
%% Use of intensity:
% The intensity of the clumped object is normalized and uniformized over
% each sub-object. The uniformization is realized by the morphological
% reconstruction of the erosion.
% 
% I = imfilter(I,fspecial('gaussian',7,1),'symmetric'); S =
% imreconstruct(imerode(I,strel('disk',EROSION_SIZE))),I,4); S =
% imopen(S,strel('disk',5));
% 
% The normalization is accomplished with dividing the image by the local
% mean
%
% S = S./imfilter(S,fspecial('gaussian',21,3),'symmetric');
%
% Finally, this resulting image is scaled for use as a multiplicative
% factor to the original confining potential (found by the distance
% transform)
%
% bounds = prctile(S(:),[10,30]); S = (S - bounds(1)) / (bounds(2) -
% bounds(1)); V = V.*S;

%% Modification of interparticle potential
% The 1/r potential has a divergence as r->0. At later times in the
% particle simulation (when the time steps are larger) if two particles get
% very close to each other than there is a very large force and the
% particle may jump out of the potential range. For this reason, and to
% increase the stability in general, the 1/r potential was modified to be
%
% 1/(r+0.2)
%
% So that the potential does not be come extreamly large as r->0.
%
% Because of this modification, the potential paramters were also resolved
% for.
