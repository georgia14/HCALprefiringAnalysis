// -----------------------------------------------------------------
void calIds(float& eta, float &phi, int& ieta, int& iphi) 
 {
   ieta = etaBin(eta);
   if(useTowerCenterEtaPhi_) eta = etaBinCenter(ieta);
   assert(ieta != 0);
   assert(ieta <= 41);
   assert(ieta >= -41);
   const static float dPhi =  2 * M_PI / 72;
   if(phi < 0) phi += 2 * M_PI;
   iphi = (int)(phi / dPhi);
   iphi++;
   if(iphi > 72) iphi = 72; 
   if(useTowerCenterEtaPhi_) {
     phi = (iphi-0.5) * dPhi;
   }
   assert(phi < 2 * M_PI);
   assert(iphi > 0);
   assert(iphi <= 72);
 }
