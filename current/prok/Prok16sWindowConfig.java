package prok;

import java.util.ArrayList;
import java.util.Locale;

/** Run-local 16S window controls and physical-edge rescue overrides.
 * @author Brian, Raiden */
public final class Prok16sWindowConfig {
	public boolean parse(String key, String value){
		if(key.equalsIgnoreCase("rrna16spartialminlen")){
			partialMinLen=CallGenes.parseSweepInt(key, value, 0, 1400);configured=true;return true;
		}
		if(key.equalsIgnoreCase("rrna16swindow")){
			if(value==null){throw new IllegalArgumentException("rrna16swindow requires legacy or containment");}
			final String v=value.toLowerCase(Locale.ROOT);
			if(!v.equals("legacy") && !v.equals("containment")){
				throw new IllegalArgumentException("rrna16swindow requires legacy or containment: "+value);
			}
			mode=v;configured=true;return true;
		}
		return false;
	}

	public void apply(ArrayList<NcrnaFamily> families){
		if(!configured){return;}
		if(!CallGenes.seededRrnaFamily(ProkObject.r16S)){
			throw new IllegalArgumentException("rrna16swindow requires rrna17=16s or rrna17=prok");
		}
		NcrnaFamily selected=null;
		for(NcrnaFamily f:families){if(f.name.equals("16S")){
			if(selected!=null){throw new IllegalArgumentException("Exactly one16S family may own window controls");}selected=f;
		}}
		if(selected==null){throw new IllegalArgumentException("Window controls have no loaded16S family");}
		final NcrnaFamily f=selected;
		final int floor=partialMinLen==null ? f.prok16sPartialMinLen : partialMinLen;
		if(floor>f.minLen){throw new IllegalArgumentException("16S partial floor cannot exceed ordinary minLen="+f.minLen);}
		if(f.kLong!=17 || f.outputType!=ProkObject.r16S || f.voteWindows || f.voteEnds || f.voteBeforePadding
				|| f.trimAlignmentExtent || !f.reuseConsensusAlignment
				|| f.boundary5NetTemplate!=null || f.boundary3NetTemplate!=null || f.rrnaEndpointFeatures!=null){
			throw new IllegalArgumentException("16S window controls require registered Quantum alignment before endpoint resources are installed");
		}
		if(mode!=null){f.preserveContainingWindows=mode.equals("containment");}
		f.prok16sPartialMinLen=floor;
		assert(!f.voteEnds && !f.trimAlignmentExtent && f.reuseConsensusAlignment)
			: "Window-only controls must not change the accepted alignment endpoints";
		System.err.println("16S window controls: preserveContaining="+f.preserveContainingWindows
			+" partialMinLen="+f.prok16sPartialMinLen);
	}
	private String mode=null;
	private boolean configured=false;
	private Integer partialMinLen=null;
}
