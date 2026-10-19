package prok;

import parse.Parse;

/** Calibrated 5.8S boundary geometry and explicit overrides for training/inference.
 * @author Raiden */
public final class R58BoundaryConfig {
	public boolean parse(String key, String value){
		if(key.equalsIgnoreCase("r58boundarystartrange")){start=range(value);return true;}
		if(key.equalsIgnoreCase("r58boundarystoprange")){stop=range(value);return true;}
		if(key.equalsIgnoreCase("r58boundaryraw")){raw=Parse.parseBoolean(value);rawExplicit=true;return true;}
		return false;
	}
	public void validate(boolean enabled){
		if(!enabled && (start!=null || stop!=null || rawExplicit && raw)){
			throw new IllegalArgumentException("Explicit 5.8S boundary ranges/raw anchors require r58boundarynn=t");
		}
	}
	public int[] offsets(boolean fivePrime, int[] defaults){
		assert(defaults!=null && defaults.length>0):"The family factory supplies the historical fallback geometry";
		final int[] chosen=fivePrime?start:stop;return chosen==null?defaults:chosen.clone();
	}
	public void apply(NcrnaFamily family, boolean enabled){
		if(family==null || !family.name.equals("r58")){throw new IllegalArgumentException("5.8S controls cannot alter another family");}
		// Vector generation loads the same ranges without loading a network.
		// CallGenes.validateNcrnaGateCombo enforces the inference-only enable gate.
		family.boundaryOnRawEndpoints=enabled && raw;
	}
	public void clear(){start=stop=null;raw=true;rawExplicit=false;}
	/** Selected intact-gene correction coverage plus two bases of margin. */
	public static int[] defaultOffsets(boolean fivePrime){
		final int low=fivePrime?-11:-46, high=fivePrime?56:37;
		final int[] offsets=new int[high-low+1];for(int i=0; i<offsets.length; i++){offsets[i]=low+i;}return offsets;
	}
	public static int[] range(String text){
		if(text==null){throw new IllegalArgumentException("Boundary range requires minimum,maximum");}
		final String[] f=text.split(",", -1);
		if(f.length!=2){throw new IllegalArgumentException("Boundary range requires exactly minimum,maximum");}
		final int low=Integer.parseInt(f[0]), high=Integer.parseInt(f[1]);
		if(low>0 || high<0 || low< -256 || high>256 || high-low>256){
			throw new IllegalArgumentException("Boundary range must include zero, stay within +/-256 and contain at most257 positions: "+text);
		}
		final int[] offsets=new int[high-low+1];for(int i=0; i<offsets.length; i++){offsets[i]=low+i;}return offsets;
	}
	private int[] start,stop;
	private boolean raw=true,rawExplicit;
}
