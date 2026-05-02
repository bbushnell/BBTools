package cardinality;

public class TestBCDLL3 {
	public static void main(String[] args){
		final int bk=60;
		final int k=31;
		final long seed=12345;
		final CompressedDynamicLogLog3 c=new CompressedDynamicLogLog3(bk, k, seed, 0);
		final BankedCompressedDynamicLogLog3 b=new BankedCompressedDynamicLogLog3(bk, k, seed, 0);

		final int cBk=c.actualBuckets(), bBk=b.actualBuckets();
		System.err.println("CDLL3 actualBuckets="+cBk+"  BCDLL3 actualBuckets="+bBk);

		final long maxCard=500000;
		boolean printed=false;

		for(long i=0; i<maxCard; i++){
			c.hashAndStore(i);
			b.hashAndStore(i);

			double[] estC=c.rawEstimates();
			double[] estB=b.rawEstimates();
			double dlcC=estC[11], dlcB=estB[11];

			if(dlcB<dlcC-0.01 && !printed){
				printed=true;
				System.err.println("\n!!! BCDLL3 DLC < CDLL3 DLC at i="+i+" trueCard="+(i+1));
				System.err.println("  DLC: C="+String.format("%.6f",dlcC)+" B="+String.format("%.6f",dlcB));
				System.err.println("  cMinZ="+c.getMinZeros()+" bMinZ="+b.getMinZeros());

				int[] nlzC=c.getNlzSnapshot();
				int[] nlzB=b.getNlzSnapshot();

				System.err.println("CDLL3 nlzCounts (full):");
				printFullNlz(nlzC);
				System.err.println("BCDLL3 nlzCounts (full):");
				printFullNlz(nlzB);

				// Check nlzCounts[0] specifically
				System.err.println("nlz[-1] (empties): C="+nlzC[0]+" B="+nlzB[0]);

				// Total counts
				int sumC=0, sumB=0;
				for(int t=0; t<66; t++){sumC+=nlzC[t]; sumB+=nlzB[t];}
				System.err.println("Total: C="+sumC+" B="+sumB+" (should both be "+cBk+")");

				// Per-bucket state
				System.err.println("Per-bucket diff:");
				for(int bi=0; bi<cBk; bi++){
					int cs=c.readBucketAt(bi);
					int bs=b.readBucketAt(bi);
					int bank=b.readBankAt(bi/10);
					int cAbs=(cs==0 && c.getMinZeros()==0) ? -1 : (cs==0 ? c.getMinZeros()-1 : (cs-1)+c.getMinZeros());
					int bAbs=(bs==0 && b.getMinZeros()+bank==0) ? -1 : (bs==0 ? b.getMinZeros()+bank-1 : (bs-1)+b.getMinZeros()+bank);
					if(cAbs!=bAbs){
						System.err.println("  bucket "+bi+": C_stored="+cs+" C_abs="+cAbs
							+"  B_stored="+bs+" B_bank="+bank+" B_abs="+bAbs);
					}
				}

				b.printDiagnostics();
				int[] bd=b.bankDistribution();
				System.err.println("Banks: b0="+bd[0]+" b1="+bd[1]+" b2="+bd[2]+" b3="+bd[3]);

				// Print all estimates for comparison
				String[] names={"Mean","HMean","HMeanM","GMean","HLL","LC",
					"Hybrid","HybDLC50","DThHyb","LCmin","DLCPure","DLC"};
				for(int e=0; e<12; e++){
					if(estC[e]!=estB[e]){
						System.err.println("  est["+e+"] "+names[e]+": C="+String.format("%.6f",estC[e])
							+" B="+String.format("%.6f",estB[e])+" diff="+String.format("%.6f",estB[e]-estC[e]));
					}
				}
				return;
			}
		}
		System.err.println("BCDLL3 DLC never fell below CDLL3 DLC in "+maxCard+" elements.");
		b.printDiagnostics();
	}

	static void printFullNlz(int[] nlz){
		StringBuilder sb=new StringBuilder("  ");
		for(int t=0; t<Math.min(nlz.length, 30); t++){
			sb.append("["+(t-1)+"]="+nlz[t]+" ");
		}
		System.err.println(sb);
	}
}
