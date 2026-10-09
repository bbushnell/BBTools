package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Locale;
import align2.MultiStateAligner9PacBio;
import align2.PacBioScoreParameters;
import dna.Data;
import map.LongHashSet;
import map.ObjectSet;
import parse.Parse;

/** Default-off, resource-bound 18S entry point. Enabling the family selects the
 * D27/D37 development recipe; explicit controls retain experimental alternatives.
 * Substitution costs stay at the existing constants until the owner selects them.
 * @author Raiden
 */
final class Euk18sRuntimeConfig {
	boolean parse(String key,String value){
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
	void apply(ArrayList<NcrnaFamily> families){
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
		final String libPath=resource(consensus,"euk18S_consensus_sequence.fa");
		final LongHashSet keys=ProkObject.loadLongKmers(resource(seeds,"euk18S_17mers.fa.gz"),ConservedRnaSeedIndex.K);
		require(keys!=null&&keys.size()>=2,"euk18S requires at least two distinct forward 17-mer keys");
		final byte[][] refs=TrnaConsensusBuilder.loadLibrary(libPath);
		final String[] names=TrnaConsensusBuilder.lastLoadedNames;
		require(refs!=null&&refs.length>0&&names!=null&&names.length==refs.length,"18S consensus must bind every model name");
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
		final NcrnaFamily f=new NcrnaFamily(FAMILY,refs,null,names,keys,17,minLen,resolvedPad,
			9,topN>0?topN:refs.length,indexFrac>0,indexFrac>0?indexMinHits:0f,indexFrac,0f,indexMinHits,0f,1f,1f,1f,
			1.01f,collapse,null,null,null,null,-1,-1,-1,-1,0f);
		f.outputType=ProkObject.r18S;f.seedMinHits=seedMinHits;f.maxLen=maxLen;
		f.seedDistinct=seedDistinct;
		f.scavengePass2=false;f.rankedModelFallback=true;f.strictIndexCutoff=true;
		f.indexScoreMargin=indexMargin;
		f.reuseConsensusAlignment=true;f.trimAlignmentExtent=false;
		f.pacBioConsensusAlignment=pacBio;
		f.pacBioCosts=costs;
		f.pacBioRolling=rolling;
		if(joinedPath!=null){f.joinedProposals=Euk18sJoinedProposals.load(joinedPath,names,joinedSupport);}
		if(useJoined && joinedPath==null){f.joinedPositions=SeedModelPositionResource.load(resource(joinedPositionPath,"euk18S_seed_positions.tsv"),keys,names,refs);}
		f.joinedSideSupport=joinedSupport;
		f.modelEndClipping=modelClip;
		f.modelClipRescue=modelClipRescue;
		f.modelClipRescueId=modelClipRescue && Float.isNaN(rescueId) ? .66f : rescueId;
		// No HBM library: detection uses the explicitly supplied identity table.
		final ArrayList<NcrnaFamily> one=new ArrayList<NcrnaFamily>(1);one.add(f);
		NcrnaModelThresholds.load(resource(cutoffs,"euk18S_model_cutoffs.tsv"),one);
		require(f.modelThresholds!=null,"18S has no accepted scalar fallback; every model needs a cutoff row");
		f.voteWindows=voteWindows;f.voteEnds=voteEnds;f.voteSlack=voteSlack;f.voteEndsMaxSd=voteMaxSd;
		f.voteBeforePadding=voteBeforePadding;
		if(voteEnds||voteWindows){
			final String expected=offsets==null && voteSeedSha80==null ? "b468f4f1faaf9dd6e56f" : voteSeedSha80;
			f.voteTable=SeedOffsetTable.load(resource(offsets,"euk18S_seed_offsets.tsv"),17,keys,expected);
		}
		if(infer){
			require(tables!=null&&nets!=null,"Uncalibrated 18S endpoint inference requires explicit table and network directories");
			final RrnaEndpointCallerFeatures.Resources loaded=RrnaEndpointResourceLoader.load(FAMILY,tables,names,refs);
			f.setRrnaEndpointInference(RrnaEndpointNetworkLoader.load(nets,loaded),null);
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
	private static String resource(String override,String fallback){
		final String chosen=override==null?fallback:override;require(chosen!=null,"Missing 18S resource path");
		final String p=new File(chosen).isFile()?chosen:Data.findPath(chosen.startsWith("?")?chosen:"?"+chosen,false);
		require(p!=null&&new File(p).isFile(),"Missing euk18S runtime resource: "+chosen);return p;
	}
	private static String path(String key,String value){require(value!=null&&!value.isEmpty(),key+" requires a path");return value;}
	private static int integer(String value,int minimum,String key){final int x=Parse.parseIntKMG(value);require(x>=minimum,key+" must be >="+minimum);return x;}
	private static float number(String value,float lo,float hi,String key){final float x=Float.parseFloat(value);require(Float.isFinite(x)&&x>=lo&&x<=hi,key+" outside ["+lo+","+hi+"]");return x;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final String FAMILY="euk18S";
	/** Explicit D48 package variant; review its manifest before building. */
	static final boolean DEFAULT_JOINED=true;
	boolean enabled=false;
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
	private boolean modelClipRescue=true;
	private float rescueId=Float.NaN;
	private String consensus,seeds,cutoffs,offsets,voteSeedSha80,tables,nets;
	private int minLen=916,maxLen=4000,flank=100,seedMinHits=2,indexMinHits=96,topN=-1,indexMargin=121,voteSlack=60;
	private int windowPad=3109;
	private float voteMaxSd=5f,collapse=.9f,indexFrac=.8f;
}
