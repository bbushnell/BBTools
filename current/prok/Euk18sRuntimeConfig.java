package prok;

import java.util.ArrayList;
import java.util.Locale;
import align2.MultiStateAligner9PacBio;
import align2.PacBioScoreParameters;
import map.ObjectSet;
import parse.Parse;

/** Default-off, resource-bound 18S entry point. Enabling the family selects the
 * D27/D37 development recipe; explicit controls retain experimental alternatives.
 * Substitution costs stay at the existing constants until the owner selects them.
 * @author Raiden
 */
public final class Euk18sRuntimeConfig {
	public boolean parse(String key,String value){
		final String k=key.toLowerCase(Locale.ROOT);
		if(k.equals("euk18s")){enabled=Parse.parseBoolean(value);return true;}
		if(!k.startsWith("euk18s")){return false;}
		if(k.equals("euk18sconsensus")){consensus=path(key,value);}
		else if(k.equals("euk18sseeds")){seeds=path(key,value);}
		else if(k.equals("euk18smodelcutoffs")){cutoffs=path(key,value);}
		else if(k.equals("euk18sjoinedproposals")){joinedPath=path(key,value);}
		else if(k.equals("euk18slivepositions")){joinedPositionPath=path(key,value);}
		else if(k.equals("euk18sjoined")){joined=Parse.parseBoolean(value);}
		else if(k.equals("euk18sjoinedsupport")){joinedSupport=integer(value,1,key);require(joinedSupport<=2,"Joined side support must be1or2");joinedSupportSet=true;}
		else if(k.equals("euk18srawends")){
			require("alignment".equalsIgnoreCase(value)||"vote".equalsIgnoreCase(value),"euk18srawends must be alignment or vote");voteEnds="vote".equalsIgnoreCase(value);
		}else if(k.equals("euk18svotewindows")){voteWindows=Parse.parseBoolean(value);}
		else if(k.equals("euk18svotebeforepadding")){voteBeforePadding=Parse.parseBoolean(value);}
		else if(k.equals("euk18svoteoffsets")){offsets=path(key,value);}
		else if(k.equals("euk18svoteseedsha80")){require(SeedOffsetTable.validSha80(value),"Expected a 20-hex seed-set sha80");voteSeedSha80=value;}
		else if(k.equals("euk18svoteslack")){voteSlack=integer(value,0,key);}
		else if(k.equals("euk18svotemaxsd")){voteMaxSd=number(value,0,Float.MAX_VALUE,key);}
		else if(k.equals("euk18sflank")){flank=integer(value,0,key);}
		else if(k.equals("euk18swindowpad")){windowPad=integer(value,0,key);}
		else if(k.equals("euk18sminlen")){minLen=integer(value,17,key);}
		else if(k.equals("euk18smaxlen")){maxLen=integer(value,17,key);}
		else if(k.equals("euk18sseedminhits")){seedMinHits=integer(value,1,key);}
		else if(k.equals("euk18sseeddistinct")){seedDistinct=Parse.parseBoolean(value);}
		else if(k.equals("euk18saligner")){
			require("quantum".equalsIgnoreCase(value) || "pacbio".equalsIgnoreCase(value),"euk18saligner must be quantum or pacbio");
			pacBio="pacbio".equalsIgnoreCase(value);
		}
		else if(k.equals("euk18smsadel4")){msaDel4=integer(value,-292,key);costsConfigured=true;}
		else if(k.equals("euk18smsaengine")){
			require("native".equalsIgnoreCase(value) || "rolling".equalsIgnoreCase(value), "euk18smsaengine must be native or rolling");
			rolling="rolling".equalsIgnoreCase(value);engineConfigured=true;
		}
		else if(k.equals("euk18smsadel5")){msaDel5=integer(value,-292,key);costsConfigured=true;}
		else if(k.equals("euk18smsains4")){msaIns4=integer(value,-292,key);costsConfigured=true;}
		else if(k.equals("euk18smsains")){msaIns=integer(value,-292,key);costsConfigured=true;}
		else if(k.equals("euk18smsadel")){msaDel=integer(value,-292,key);costsConfigured=true;}
		else if(k.equals("euk18smsasubscale")){msaSubScale=number(value,Float.MIN_VALUE,1,key);costsConfigured=true;}
		else if(k.equals("euk18smodelclip")){
			modelClipRescue="rescue".equalsIgnoreCase(value);
			modelClip=!modelClipRescue && Parse.parseBoolean(value);
		}
		else if(k.equals("euk18srescueid")){rescueId=number(value,0,1,key);}
		else if(k.equals("euk18sindexminhits")){indexMinHits=integer(value,1,key);}
		else if(k.equals("euk18stopn")){topN=integer(value,1,key);}
		else if(k.equals("euk18sindexmargin")){indexMargin=integer(value,-1,key);}
		else if(k.equals("euk18sindexfrac")){indexFrac=number(value,0,1,key);}
		else if(k.equals("euk18scollapsefrac")){collapse=number(value,0,1,key);}
		else if(k.equals("euk18sendpoint")){infer=Parse.parseBoolean(value);}
		else if(k.equals("euk18sendpointtables")){tables=path(key,value);}
		else if(k.equals("euk18sendpointnets")){nets=path(key,value);}
		else{return false;}
		configured=true;return true;
	}

	/** Validate and construct privately; no partial family is published on failure. */
	public void apply(ArrayList<NcrnaFamily> families){
		require(families!=null,"18S registration requires the caller family list");
		if(!enabled){require(!configured,"euk18S resource/control options require euk18s=t");return;}
		for(NcrnaFamily f:families){require(!FAMILY.equals(f.name)&&!"18S".equals(f.name),"Duplicate 18S registration; choose euk18s=t or the historical rrna17 profile");}
		require(maxLen>=minLen,"euk18smaxlen must be at least euk18sminlen");
		require(!pacBio || (!modelClip && !infer && !voteEnds),"Experimental PacBio primary requires no primary clipping, endpoint NN off and raw alignment ends");
		require(!costsConfigured || pacBio,"Experimental MSA costs require euk18saligner=pacbio");
		require(!engineConfigured || pacBio, "euk18smsaengine requires euk18saligner=pacbio");
		require(joinedPath==null || joinedPositionPath==null,"Choose saved proposal replay or live discovery, never both");
		final boolean suppliedJoined=joinedPath!=null || joinedPositionPath!=null;
		final boolean useJoined=joined!=null ? joined.booleanValue() : suppliedJoined || (DEFAULT_JOINED && pacBio && rolling);
		require(useJoined || !suppliedJoined,"euk18sjoined=f conflicts with an explicit joined proposal resource");
		require(!joinedSupportSet || useJoined,"euk18sjoinedsupport requires joined discovery or replay");
		require(!useJoined || (pacBio && rolling && !infer && !voteEnds && minLen==916),"Joined route requires rolling PacBio, raw alignment ends, no endpoint NN and916core minimum; use euk18sjoined=f for Quantum");
		final PacBioScoreParameters costs=PacBioScoreParameters.experimental(msaDel4,msaDel5,msaIns4,msaIns,msaDel,msaSubScale);
		require(Float.isNaN(rescueId) || modelClipRescue,"euk18srescueid requires euk18smodelclip=rescue");
		require(!(modelClip || modelClipRescue) || (!infer && !voteEnds),"Experimental model clipping requires endpoint NN off and raw alignment ends");
		require(!voteBeforePadding || (voteWindows && windowPad>=0 && collapse>0),"Vote-before-padding requires window votes, an explicit measured windowpad, and positive core overlap");
		require(indexMargin<0 || topN<0,"Choose euk18sindexmargin or euk18stopn; a score-margin shortlist has no topN cap");
		require(indexFrac==0 || topN<0,"euk18sindexfrac uses all models; do not specify euk18stopn");
		require(infer||(tables==null&&nets==null),"Endpoint directories require euk18sendpoint=t");
		require((voteEnds||voteWindows)||(offsets==null&&voteSeedSha80==null),"Vote offsets require voted windows or raw ends");
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.hbm=false;options.endpoints=false;options.cutoffs=false;
		options.resources.consensus=consensus;options.resources.seeds=seeds;options.resources.cutoffs=cutoffs;
		options.resources.voteOffsets=offsets;options.resources.joinedModelPositions=joinedPositionPath;
		options.resources.tables=tables;options.resources.endpoints=nets;
		options.minLen=minLen;options.windowPad=windowPad;options.seedMinHits=seedMinHits;
		if(windowPad<0){options.longestPadFlank=flank;}
		options.indexK=9;options.topN=topN;options.minHits=indexMinHits;
		options.adaptive=indexFrac>0;options.floor=indexFrac>0 ? indexMinHits : 0f;
		options.topFrac=indexFrac;options.qFrac=0f;options.collapse=collapse;
		final NcrnaFamily f=NcrnaFamilyLoader.load(FAMILY,options,resolver);
		final byte[][] refs=f.library;
		final String[] names=f.modelNames;
		require(refs!=null&&refs.length>0&&names!=null&&names.length==refs.length,"18S consensus must bind every model name");
		if(consensus==null){validateDefaultNames(names);}
		final ObjectSet<String> unique=new ObjectSet<String>(String.class);int longest=0;
		for(int i=0;i<refs.length;i++){
			require(names[i]!=null&&names[i].matches("[A-Za-z0-9_]+")&&unique.add(names[i]),"18S model names must be unique safe tokens");
			require(refs[i]!=null&&refs[i].length>=17,"18S model is shorter than the shared seed K");longest=Math.max(longest,refs[i].length);
		}
		// Existing window builder starts at seedCentre-pad and ends at
		// seedCentre+pad+17. A full longest-model span plus flank on each side
		// of every seed avoids assuming where within a ~1.8 kb gene that seed lies.
		require(longest<=Integer.MAX_VALUE-flank-17,"18S window length would overflow integer coordinates");
		final int resolvedPad=windowPad<0?longest+flank:windowPad;
		require(resolvedPad<=Integer.MAX_VALUE-17,"18S fallback window padding must leave room for its seed span");
		require(f.windowPad==resolvedPad,"Shared loader must preserve longest-model plus flank or explicit 18S padding");
		f.outputType=ProkObject.r18S;f.seedMinHits=seedMinHits;f.maxLen=maxLen;
		f.seedDistinct=seedDistinct;
		f.scavengePass2=false;f.rankedModelFallback=true;f.strictIndexCutoff=true;
		f.indexScoreMargin=indexMargin;
		f.reuseConsensusAlignment=true;f.trimAlignmentExtent=false;
		f.pacBioConsensusAlignment=pacBio;
		f.pacBioCosts=costs;
		f.pacBioRolling=rolling;
		if(joinedPath!=null){f.joinedProposals=Euk18sJoinedProposals.load(joinedPath,names,joinedSupport);}
		if(useJoined && joinedPath==null){NcrnaFamilyLoader.loadJoinedPositions(f,options.resources,resolver);}
		f.joinedSideSupport=joinedSupport;
		f.modelEndClipping=modelClip;
		f.modelClipRescue=modelClipRescue;
		f.modelClipRescueId=modelClipRescue && Float.isNaN(rescueId) ? .66f : rescueId;
		// No HBM library: detection uses the explicitly supplied identity table.
		NcrnaFamilyLoader.loadCutoffs(f,options.resources,resolver);
		require(f.modelThresholds!=null,"18S has no accepted scalar fallback; every model needs a cutoff row");
		f.voteWindows=voteWindows;f.voteEnds=voteEnds;f.voteSlack=voteSlack;f.voteEndsMaxSd=voteMaxSd;
		f.voteBeforePadding=voteBeforePadding;
		if(voteEnds||voteWindows){
			final String expected=offsets==null && voteSeedSha80==null ? "47bd7e395c8398204083" : voteSeedSha80;
			NcrnaFamilyLoader.loadVoteOffsets(f,options.resources,resolver,expected);
		}
		if(infer){
			require(tables!=null&&nets!=null,"Uncalibrated 18S endpoint inference requires explicit table and network directories");
			f.setRrnaEndpointInference(NcrnaFamilyLoader.loadEndpoints(f,options.resources,resolver),null);
		}
		families.add(f);
		System.err.println("euk18S controls: models="+refs.length+" seedK=17 seedMinHits="+seedMinHits
			+" seedMode="+(seedDistinct?"distinct":"occurrence")
			+" indexK=9 topN="+(indexMargin<0 && indexFrac==0 ? Integer.toString(f.indexTopN) : "uncapped")
			+" indexMargin="+indexMargin+" indexMinHits="+indexMinHits+" windowPad="+f.windowPad
			+" flank="+flank+" minLen="+minLen+" maxLen="+maxLen+" rawEnds="+(voteEnds?"vote":"alignment")
			+" voteWindows="+voteWindows+" voteBeforePadding="+voteBeforePadding+" voteSlack="+voteSlack+" voteMaxSd="+voteMaxSd+" endpoint="+infer+" indexFrac="+indexFrac+" modelClip="+modelClip+" modelClipRescue="+modelClipRescue
			+" rescueId="+(Float.isNaN(f.modelClipRescueId)?"model-cutoff":Float.toString(f.modelClipRescueId))+" aligner="+(pacBio?"pacbio":"quantum")
			+" msaDel4="+costs.del4+" msaDel5="+costs.del5+" msaIns4="+costs.ins4+" msaIns="+costs.ins+" msaDel="+costs.del
			+" msaSubScale="+msaSubScale+" msaSub="+costs.sub+" msaSubR="+costs.subR+" msaSub2="+costs.sub2+" msaSub3="+costs.sub3
			+" msaEngine="+(rolling ? "rolling" : "native")+" joined="+useJoined+" joinedSupport="+joinedSupport);
	}
	static void validateDefaultNames(String[] names){
		require(names!=null && names.length==DEFAULT_MODEL_NAMES.length, "Default euk18S library requires all ten declared models");
		for(int i=0; i<names.length; i++){
			require(DEFAULT_MODEL_NAMES[i].equals(names[i]), "Default euk18S model roster mismatch at "+i+": expected="+DEFAULT_MODEL_NAMES[i]+" observed="+names[i]);
		}
	}
	private static String path(String key,String value){require(value!=null&&!value.isEmpty(),key+" requires a path");return value;}
	private static int integer(String value,int minimum,String key){final int x=Parse.parseIntKMG(value);require(x>=minimum,key+" must be >="+minimum);return x;}
	private static float number(String value,float lo,float hi,String key){final float x=Float.parseFloat(value);require(Float.isFinite(x)&&x>=lo&&x<=hi,key+" outside ["+lo+","+hi+"]");return x;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String FAMILY="euk18S";
	private static final String[] DEFAULT_MODEL_NAMES={"euk18S_universal", "euk18S_fungi", "euk18S_plant",
		"euk18S_invertebrate", "euk18S_vertebrate_other", "euk18S_vertebrate_mammalian",
		"euk18S_mitochondria", "euk18S_protist", "euk18S_protozoa_component_001", "euk18S_protozoa_component_003"};
	/** Explicit D48 package variant; review its manifest before building. */
	static final boolean DEFAULT_JOINED=true;
	public boolean enabled=false;
	private boolean configured=false,voteWindows=true,voteEnds=false,infer=false,voteBeforePadding=true;
	private boolean seedDistinct=true;
	private boolean modelClip=false;
	private boolean pacBio=true;
	private boolean costsConfigured=false;
	private boolean rolling=true, engineConfigured=false;
	private String joinedPath,joinedPositionPath;
	private Boolean joined=null;
	private int joinedSupport=1;
	private boolean joinedSupportSet=false;
	private int msaDel4=MultiStateAligner9PacBio.POINTS_DEL4,msaDel5=MultiStateAligner9PacBio.POINTS_DEL5,
		msaIns4=MultiStateAligner9PacBio.POINTS_INS4,msaIns=MultiStateAligner9PacBio.POINTS_INS,msaDel=MultiStateAligner9PacBio.POINTS_DEL;
	private float msaSubScale=.66f;
	private boolean modelClipRescue=false;
	private float rescueId=Float.NaN;
	private String consensus,seeds,cutoffs,offsets,voteSeedSha80,tables,nets;
	private int minLen=916,maxLen=4000,flank=100,seedMinHits=2,indexMinHits=96,topN=-1,indexMargin=121,voteSlack=60;
	private int windowPad=3109;
	private float voteMaxSd=5f,collapse=.9f,indexFrac=.8f;
}
