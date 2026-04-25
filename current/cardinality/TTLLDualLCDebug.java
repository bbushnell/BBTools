package cardinality;

import shared.Tools;
import java.lang.reflect.Field;

/**
 * Diagnostic: dump TTLL register state and compare LC vs DualLC at each cardinality.
 * Usage: java cardinality.TTLLDualLCDebug [buckets=16] [maxcard=60]
 */
public class TTLLDualLCDebug {

	public static void main(String[] args) throws Exception {
		int buckets=16, maxCard=60;
		for(String arg : args){
			String[] ab=arg.split("=");
			if(ab[0].equals("buckets")){buckets=Integer.parseInt(ab[1]);}
			else if(ab[0].equals("maxcard")){maxCard=Integer.parseInt(ab[1]);}
		}

		TwinTailLogLog.ENCODING_MODE=0;
		TwinTailLogLog ttll=new TwinTailLogLog(buckets, 31, 12345, 0);

		// Reflective access to private fields
		Field regsF=TwinTailLogLog.class.getDeclaredField("regs");
		regsF.setAccessible(true);
		Field geF=TwinTailLogLog.class.getDeclaredField("globalExp");
		geF.setAccessible(true);

		System.out.println("Card\tLC\tDualLC\tB\tV_lc\tVirtB\tVirtV\tGE\tRegisters");

		for(int card=1; card<=maxCard; card++){
			ttll.add(card);

			byte[] regs=(byte[])regsF.get(ttll);
			int globalExp=geF.getInt(ttll);

			// Standard LC: from register exponents only
			int filled=0;
			StringBuilder regStr=new StringBuilder();
			for(int i=0; i<buckets; i++){
				int b=regs[i]&0xFF;
				int le=(b>>>4)&0xF;
				int absNlz=globalExp+le;
				int h0=b&0x3;
				int h1=(b>>>2)&0x3;
				if(b!=0 || globalExp>0) filled++;
				if(i>0) regStr.append(' ');
				regStr.append(absNlz).append('.');
				regStr.append(Integer.toBinaryString(4|h1).substring(1));
				regStr.append(Integer.toBinaryString(4|h0).substring(1));
			}
			int V_lc=buckets-filled;
			double lcEst=(V_lc>0) ? (double)buckets*Math.log((double)buckets/V_lc) : 0;

			// DualLC
			double dualEst=ttll.dualBucketLC();

			// Recompute virtual counts for display
			int virtualTotal=0, virtualEmpty=0;
			for(int i=0; i<buckets; i++){
				int b=regs[i]&0xFF;
				if(b==0 && globalExp==0){
					virtualTotal+=2; virtualEmpty+=2; continue;
				}
				int le=(b>>>4)&0xF;
				int absNlz=globalExp+le;
				for(int t=0; t<2; t++){
					int tail=(t==0) ? (b&0x3) : ((b>>>2)&0x3);
					if(tail>=2){ virtualTotal++; }
					else if(tail==1){ virtualTotal++; }
					else{
						if(absNlz<2){ virtualTotal++; virtualEmpty++; }
					}
				}
			}

			System.out.printf("%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%s%n",
				card, lcEst, dualEst, buckets, V_lc, virtualTotal, virtualEmpty, globalExp, regStr);
		}
	}
}
