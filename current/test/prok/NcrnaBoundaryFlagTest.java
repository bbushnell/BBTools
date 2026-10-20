package test.prok;

import java.util.Arrays;
import prok.CallGenes;
import prok.GeneCaller;
import prok.NcrnaFamily;

/** Verifies CLI-to-loaded-network selection in a fresh JVM with every prok family retained.
 * @author Ganyu */
public final class NcrnaBoundaryFlagTest {

	public static void main(String[] args){
		if(args.length<3 || !args[0].startsWith("expectrnasep=") || !args[1].startsWith("expecttmrna=")){
			throw new IllegalArgumentException("expectrnasep=t|f expecttmrna=t|f followed by ordinary CallGenes arguments");
		}
		final boolean rnasep=expected(args[0]), tmrna=expected(args[1]);
		boolean small=true, large=true, rf13=true, rf85=true;int mask=63, first=2;
		int[] start13={-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8}, stop13=start13;
		int[] start85={-1,0,1}, stop85=start85;
		while(first<args.length && args[first].startsWith("expect")){
			final String arg=args[first++];
			if(arg.startsWith("expectsmall=")){small=expected(arg);}
			else if(arg.startsWith("expectlarge=")){large=expected(arg);}
			else if(arg.startsWith("expectmask=")){mask=Integer.parseInt(arg.substring(11));}
			else if(arg.startsWith("expect13=")){rf13=expected(arg);}
			else if(arg.startsWith("expect85=")){rf85=expected(arg);}
			else if(arg.startsWith("expect13start=")){start13=expectedOffsets(arg);}
			else if(arg.startsWith("expect13stop=")){stop13=expectedOffsets(arg);}
			else if(arg.startsWith("expect85start=")){start85=expectedOffsets(arg);}
			else if(arg.startsWith("expect85stop=")){stop85=expectedOffsets(arg);}
			else{throw new IllegalArgumentException("Unknown fixture expectation: "+arg);}
		}
		require(mask>=0 && mask<=63, "Expected generic-family mask must use only the six declared bits");
		new CallGenes(Arrays.copyOfRange(args, first, args.length));
		int seen=0;
		for(NcrnaFamily family:GeneCaller.ncrnaFamilies){
			final boolean selected;
			switch(family.name){
				case "rnasep":selected=rnasep;seen|=1;break;
				case "tmrna":selected=tmrna;seen|=2;break;
				case "srp_small":selected=small;seen|=4;break;
				case "srp_large":selected=large;seen|=8;break;
				case "sixs_rf00013":selected=rf13;seen|=16;offsets(family, start13, stop13);break;
				case "sixs_rf01685":selected=rf85;seen|=32;offsets(family, start85, stop85);break;
				default:continue;
			}
			require((family.boundary5NetTemplate!=null)==selected
				&& (family.boundary3NetTemplate!=null)==selected
				&& (family.boundaryStartTable!=null)==selected
				&& (family.boundaryStopTable!=null)==selected,
				"CLI selection did not reach complete boundary resources for "+family.name);
		}
		require(seen==mask, "Boundary defaults must retain exactly the selected caller families: expected="+mask+" observed="+seen);
		// Explicit false must also defeat a true global value; unrelated families inherit it.
		CallGenes.SIXS_RF00013_BOUNDARY_CONFIG.clear();CallGenes.SIXS_RF01685_BOUNDARY_CONFIG.clear();
		CallGenes.NCRNA_BOUNDARY_NN_ENABLED=true;
		CallGenes.RNASEP_BOUNDARY_NN_OVERRIDE=false;CallGenes.TMRNA_BOUNDARY_NN_OVERRIDE=false;
		CallGenes.SRPSMALL_BOUNDARY_NN_OVERRIDE=null;CallGenes.SRPLARGE_BOUNDARY_NN_OVERRIDE=null;
		require(!CallGenes.ncrnaBoundaryEnabled("rnasep") && !CallGenes.ncrnaBoundaryEnabled("tmrna")
			&& CallGenes.ncrnaBoundaryEnabled("srp_small") && CallGenes.ncrnaBoundaryEnabled("sixs_rf00013"),
			"Explicit opt-outs must not change global inheritance for unrelated families");
		CallGenes.SRPSMALL_BOUNDARY_NN_OVERRIDE=false;CallGenes.SRPLARGE_BOUNDARY_NN_OVERRIDE=false;
		require(!CallGenes.ncrnaBoundaryEnabled("srp_small") && !CallGenes.ncrnaBoundaryEnabled("srp_large")
			&& CallGenes.ncrnaBoundaryEnabled("sixs_rf00013"), "SRP opt-outs must leave global6S strict behavior unchanged");
		System.out.println("BOUNDARY_FLAG_PASS rnasep="+rnasep+" tmrna="+tmrna+" srp_small="+small+" srp_large="+large+" generic_mask="+seen);
		System.out.println("SIXS_BOUNDARY_FLAG_PASS rf00013="+rf13+" rf01685="+rf85+" loaded_ranges=true");
	}
	private static void offsets(NcrnaFamily family, int[] start, int[] stop){
		require(Arrays.equals(start, family.boundaryStartOffsets) && Arrays.equals(stop, family.boundaryStopOffsets),
			"CLI correction ranges must reach the loaded scorer input for "+family.name);
	}
	private static int[] expectedOffsets(String arg){
		final String[] fields=arg.substring(arg.indexOf('=')+1).split(",", -1);
		final int[] result=new int[fields.length];
		for(int i=0; i<fields.length; i++){result[i]=Integer.parseInt(fields[i]);}
		require(result.length>0, "Fixture must name every expected offset explicitly");return result;
	}
	private static boolean expected(String arg){
		final String value=arg.substring(arg.indexOf('=')+1);
		if("t".equals(value)){return true;}
		if("f".equals(value)){return false;}
		throw new IllegalArgumentException("Expected fixture values must be t or f: "+arg);
	}
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
