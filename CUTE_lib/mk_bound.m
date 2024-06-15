function[positions]=mk_bound(x1,x2,z1,z2)
 positions=[x1,0,z1;x1,0,z2;x2,0,z1;x2,0,z2];
 xmid=(x1+x2)/2;
 zmid=(z1+z2)/2;
 positions=cat(1,positions,[xmid,0,z1;xmid,0,z2;x1,0,zmid;x2,0,zmid]);
end