// phi is already relative; settling is a relative slip velocity.
// Rebuild before transport; never subtract meshPhi twice from stored phip.
phip = phi + fvc::flux(Ws);
