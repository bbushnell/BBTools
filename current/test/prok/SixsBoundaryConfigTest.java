package test.prok;

import java.util.Arrays;
import prok.NcrnaFamilyLoader;
import prok.SixsBoundaryConfig;

/** Exercises both subtype configurations without I/O, global state or caller-policy changes.
 * @author Ganyu */
public final class SixsBoundaryConfigTest {

	public static void main(String[] args){
		for(String family:new String[]{"sixs_rf00013","sixs_rf01685"}){exercise(family);}
		fails(()->new SixsBoundaryConfig("sixs"));
		fails(()->new SixsBoundaryConfig("rnasep"));
		fails(()->new SixsBoundaryConfig(null));
		System.out.println("SIXS_BOUNDARY_CONFIG_TEST_PASS subtypes=2 isolation=true order=true gates=true defensive_copies=true");
	}

	private static void exercise(String family){
		final String prefix=family.replace("_","");
		final String other=family.equals("sixs_rf00013")?"sixs_rf01685":"sixs_rf00013";
		final String otherPrefix=other.replace("_","");
		final SixsBoundaryConfig c=new SixsBoundaryConfig(family),sibling=new SixsBoundaryConfig(other);
		check(!c.enabled(false) && c.enabled(true), "Absent override must inherit the existing global boolean");
		for(String suffix:SUFFIXES){
			check(!c.parse(otherPrefix+suffix,"unused"), "Another subtype's option must not be consumed");
			check(!sibling.parse(prefix+suffix,"unused"), "Subtype isolation must work in both directions");
		}
		check(!c.parse(prefix+"boundaryraw","t"), "Unknown controls must remain unconsumed");
		check(!c.enabled(false) && c.enabled(true), "Unrecognized keys cannot change enable state");

		final NcrnaFamilyLoader.Options untouched=options();
		final int[] originalStart=untouched.boundaryStartOffsets,originalStop=untouched.boundaryStopOffsets;
		c.apply(family,untouched,false);
		check(!untouched.boundary && untouched.boundaryStartOffsets==originalStart && untouched.boundaryStopOffsets==originalStop
			&& untouched.resources.boundary.equals("existing.net"), "Absent selections must preserve caller-provided resources and ranges");
		unchanged(untouched);
		c.apply(family,untouched,true);check(untouched.boundary, "Global enable must reach loader options");

		check(c.parse(prefix+"boundarynn","t"), "Explicit subtype enable must parse");
		fails(()->c.validate(false,true,false));fails(()->c.validate(true,false,false));fails(()->c.validate(false,false,true));
		c.validate(true,true,false);check(c.enabled(false), "Explicit true overrides a false global flag");
		c.parse(prefix+"boundarynn","f");
		c.validate(false,false,true);check(!c.enabled(true), "Explicit false opts out of a true global flag");
		final NcrnaFamilyLoader.Options off=options();off.boundary=true;c.apply(family,off,true);
		check(!off.boundary, "Explicit false must reach loader options");unchanged(off);

		c.clear();
		for(String suffix:new String[]{"boundarynet","boundarystarttable","boundarystoptable","boundarystartrange","boundarystoprange"}){
			final String value=suffix.endsWith("range")?"-1,1":"not-opened.data";
			c.clear();check(c.parse(prefix+suffix,value), "Every explicit resource/range control must parse");
			fails(()->c.validate(true,true,false));
			final NcrnaFamilyLoader.Options rejected=options();
			fails(()->c.apply(family,rejected,false));
			check(!rejected.boundary && rejected.resources.boundary.equals("existing.net"),
				"Rejected application must not partially mutate loader state");
			fails(()->c.validate(false,true,true));
			fails(()->c.validate(true,false,true));
			fails(()->c.validate(false,false,true));
			c.validate(true,true,true);
			c.parse(prefix+"boundarynn","f");fails(()->c.validate(true,true,true));
			c.parse(prefix+"boundarynn","t");c.validate(true,true,false);
		}

		final SixsBoundaryConfig early=new SixsBoundaryConfig(family),late=new SixsBoundaryConfig(family);
		configure(early,prefix,true);configure(late,prefix,false);
		early.validate(true,true,false);late.validate(true,true,false);
		final NcrnaFamilyLoader.Options a=options(),b=options();
		early.apply(family,a,false);late.apply(family,b,false);
		check(a.boundary && b.boundary && a.resources.boundary.equals("new.net")
			&& a.resources.boundaryStart.equals("start.table") && a.resources.boundaryStop.equals("stop.table"),
			"Explicit paths must reach only their designated loader roles");
		check(a.resources.boundary.equals(b.resources.boundary) && a.resources.boundaryStart.equals(b.resources.boundaryStart)
			&& a.resources.boundaryStop.equals(b.resources.boundaryStop)
			&& Arrays.equals(a.boundaryStartOffsets,b.boundaryStartOffsets) && Arrays.equals(a.boundaryStopOffsets,b.boundaryStopOffsets),
			"Enable-before-resources and enable-after-resources must yield identical selections");
		check(Arrays.equals(a.boundaryStartOffsets,new int[]{-2,-1,0,1,2,3})
			&& Arrays.equals(a.boundaryStopOffsets,new int[]{0}), "Ranges are ordered, inclusive and zero-containing");
		unchanged(a);unchanged(b);
		a.boundaryStartOffsets[0]=99;a.boundaryStopOffsets[0]=99;
		final NcrnaFamilyLoader.Options copy=options();early.apply(family,copy,false);
		check(copy.boundaryStartOffsets[0]==-2 && copy.boundaryStopOffsets[0]==0
			&& b.boundaryStartOffsets[0]==-2 && b.boundaryStopOffsets[0]==0,
			"Each application must clone offsets rather than expose or share configuration arrays");
		fails(()->early.apply(other,options(),false));fails(()->early.apply(null,options(),false));
		fails(()->early.apply(family,null,false));

		for(String bad:new String[]{null,"","1,4","-5,-1","-300,2","-200,200","0","0,1,2","2,-2","x,1"}){
			fails(()->c.parse(prefix+"boundarystartrange",bad));
			fails(()->c.parse(prefix+"boundarystoprange",bad));
		}
		for(String suffix:new String[]{"boundarynet","boundarystarttable","boundarystoptable"}){
			fails(()->c.parse(prefix+suffix,null));fails(()->c.parse(prefix+suffix,""));fails(()->c.parse(prefix+suffix,"  "));
		}
		early.clear();early.validate(false,false,false);
		check(!early.enabled(false) && early.enabled(true), "Clear must restore nullable global inheritance");
		final NcrnaFamilyLoader.Options cleared=options();early.apply(family,cleared,false);
		check(!cleared.boundary && cleared.resources.boundary.equals("existing.net")
			&& cleared.resources.boundaryStart.equals("existing.start") && cleared.resources.boundaryStop.equals("existing.stop")
			&& Arrays.equals(cleared.boundaryStartOffsets,new int[]{-1,0}) && Arrays.equals(cleared.boundaryStopOffsets,new int[]{0,1}),
			"Clear must forget paths and ranges as well as explicit enable");
		unchanged(cleared);
	}

	private static void configure(SixsBoundaryConfig c, String prefix, boolean first){
		if(first){c.parse(prefix+"boundarynn","t");}
		c.parse(prefix+"boundarynet","new.net");c.parse(prefix+"boundarystarttable","start.table");c.parse(prefix+"boundarystoptable","stop.table");
		c.parse(prefix.toUpperCase(java.util.Locale.ROOT)+"BOUNDARYSTARTRANGE","-2,3");c.parse(prefix+"boundarystoprange","0,0");
		if(!first){c.parse(prefix+"boundarynn","t");}
	}
	private static NcrnaFamilyLoader.Options options(){
		final NcrnaFamilyLoader.Options o=new NcrnaFamilyLoader.Options();
		o.enabled=false;o.hbm=false;o.seeds=false;o.endpoints=false;o.cutoffs=false;
		o.adaptive=true;o.indexK=7;o.topN=60;o.minHits=2;o.seedMinHits=3;o.windowPad=125;o.minLen=60;o.seedK=17;o.longestPadFlank=21;
		o.floor=4;o.topFrac=.25f;o.qFrac=.1f;o.weight=2;o.idPass=.7f;o.idBorderline=.6f;o.hbmPass=.64f;
		o.scoreA=3;o.scoreB=20;o.collapse=.85f;o.boundaryMarginStart=.01f;o.boundaryMarginStop=.02f;
		o.boundaryStartOffsets=new int[]{-1,0};o.boundaryStopOffsets=new int[]{0,1};
		o.resources.root="root";o.resources.consensus="consensus";o.resources.models="models";o.resources.seeds="seeds";
		o.resources.tables="tables";o.resources.endpoints="endpoints";o.resources.cutoffs="cutoffs";
		o.resources.boundary="existing.net";o.resources.boundaryStart="existing.start";o.resources.boundaryStop="existing.stop";
		o.resources.boundary3="third";o.resources.voteOffsets="votes";o.resources.joinedModelPositions="positions";
		return o;
	}
	private static void unchanged(NcrnaFamilyLoader.Options o){
		check(!o.enabled && !o.hbm && !o.seeds && !o.endpoints && !o.cutoffs && Boolean.TRUE.equals(o.adaptive)
			&& o.indexK==7 && o.topN==60 && o.minHits==2 && o.seedMinHits==3 && o.windowPad==125 && o.minLen==60
			&& o.seedK==17 && o.longestPadFlank==21 && o.floor==4 && o.topFrac==.25f && o.qFrac==.1f && o.weight==2
			&& o.idPass==.7f && o.idBorderline==.6f && o.hbmPass==.64f && o.scoreA==3 && o.scoreB==20 && o.collapse==.85f
			&& o.boundaryMarginStart==.01f && o.boundaryMarginStop==.02f,
			"Boundary overrides must not change enable, mapping, verification, scoring, geometry or margins");
		check(o.resources.root.equals("root") && o.resources.consensus.equals("consensus") && o.resources.models.equals("models")
			&& o.resources.seeds.equals("seeds") && o.resources.tables.equals("tables") && o.resources.endpoints.equals("endpoints")
			&& o.resources.cutoffs.equals("cutoffs") && o.resources.boundary3.equals("third")
			&& o.resources.voteOffsets.equals("votes") && o.resources.joinedModelPositions.equals("positions"),
			"Boundary overrides must preserve every unrelated resource selection");
	}
	private static void fails(Runnable action){
		try{action.run();}catch(IllegalArgumentException expected){return;}
		throw new AssertionError("Invalid 6S boundary configuration was accepted");
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] SUFFIXES={"boundarynn","boundarynet","boundarystarttable","boundarystoptable","boundarystartrange","boundarystoprange"};
}
