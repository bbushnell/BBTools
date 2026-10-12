package prok;

import java.util.ArrayList;
import java.util.Locale;

/** Explicit23S-only window experiment; absent flags retain the measured collapse policy.
 * @author Brian, Raiden */
final class Prok23sWindowConfig {
	boolean parse(String key, String value){
		if(key.equalsIgnoreCase("rrna23spartialminlen")){
			partialMinLen=CallGenes.parseSweepInt(key,value,0,2600);return true;
		}
		if(!key.equalsIgnoreCase("rrna23swindow")){return false;}
		if(value==null){throw new IllegalArgumentException("rrna23swindow requires legacy or containment");}
		final String v=value.toLowerCase(Locale.ROOT);
		if(!v.equals("legacy") && !v.equals("containment")){
			throw new IllegalArgumentException("rrna23swindow requires legacy or containment: "+value);
		}
		mode=v;return true;
	}
	void apply(ArrayList<NcrnaFamily> families){
		if(mode==null && partialMinLen==null){return;}
		if(!CallGenes.RRNA17_ENABLED || !("23s".equals(CallGenes.RRNA17_PROFILE) || "prok".equals(CallGenes.RRNA17_PROFILE))){
			throw new IllegalArgumentException("rrna23swindow requires rrna17=23s or rrna17=prok");
		}
		NcrnaFamily selected=null;
		for(NcrnaFamily f:families){if(f.name.equals("23S")){
			if(selected!=null){throw new IllegalArgumentException("Exactly one23S family may own window controls");}selected=f;
		}}
		if(selected==null){throw new IllegalArgumentException("Window controls have no loaded23S family");}
		final NcrnaFamily f=selected;
		final int floor=partialMinLen==null?f.prok16sPartialMinLen:partialMinLen;
		if(floor>f.minLen){throw new IllegalArgumentException("23S partial floor cannot exceed ordinary minLen="+f.minLen);}
		if(f.kLong!=17 || f.outputType!=ProkObject.r23S || f.voteWindows || f.voteEnds || f.voteBeforePadding
				|| f.trimAlignmentExtent || !f.reuseConsensusAlignment
				|| f.boundary5NetTemplate!=null || f.boundary3NetTemplate!=null || f.rrnaEndpointFeatures!=null){
			throw new IllegalArgumentException("23S window controls require raw Quantum alignment before endpoint resources are installed");
		}
		if(mode!=null){f.preserveContainingWindows=mode.equals("containment");}
		f.prok16sPartialMinLen=floor;
		assert(!f.voteEnds && !f.trimAlignmentExtent && f.reuseConsensusAlignment)
			: "Window-only controls must preserve the accepted alignment endpoint contract";
		System.err.println("23S window controls: preserveContaining="+f.preserveContainingWindows+" partialMinLen="+floor);
	}
	private String mode=null;
	private Integer partialMinLen=null;
}
