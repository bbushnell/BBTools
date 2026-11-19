package prok;

import java.util.ArrayList;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import map.LongHashSet;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.Read;

/** Builds the fixed RF00174 member-voted profile resource once, before callers.
 * Forked from template/A_SampleByteFile's argument/validation/process lifecycle;
 * existing native member/FASTA owners handle the joined construction inputs.
 * Voting is one owner pass with reused sketch scratch, not per-caller work.
 * No genome, alignment, verification or caller invocation.
 * @author Brian, Raiden
 */
public final class BuildNcrnaDdlProfiles{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, BuildNcrnaDdlProfiles.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Config c=parse(pp.args);
			DynamicDemiLog.setExponent(4);
			build(c);
			timer.stop();
			pp.outstream.println("Profile models: "+c.models+"; kept member votes: "+c.members);
			pp.outstream.println("Time: "+timer);
			pp.outstream.println("NCRNA_DDL_PROFILE_RESOURCE_BUILD_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------    Initialization Helpers    ----------------*/
	/*--------------------------------------------------------------*/

	private static Config parse(final String[] args){
		assert(args!=null) : "PreParser must supply the expanded native argument list";
		final Config c=new Config();
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("ref")){c.references=b;}
			else if(a.equals("ledger")){c.ledger=b;}
			else if(a.equals("assignments")){c.assignments=b;}
			else if(a.equals("census")){c.census=b;}
			else if(a.equals("summary")){c.summary=b;}
			else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
			else if(a.equals("expectedmembers")){c.members=Parse.parseIntKMG(b);}
			else if(a.equals("dual")){c.dual=Parse.parseBoolean(b);}
			else if(a.equals("partitionbins")){c.partitionBins=Parse.parseIntKMG(b);}
			else if(a.equals("selected")){c.selected=Parse.parseBoolean(b);}
			else if(a.equals("expectedinput")){c.inputMembers=Parse.parseIntKMG(b);}
			else if(a.equals("paddingblacklist")){c.paddingBlacklist=b;}
			else if(a.equals("paddingfit")){c.paddingFit=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown profile-build parameter: "+arg);}
		}
		require(parser.in1!=null && parser.out1!=null && parser.in2==null && parser.out2==null
			&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
			"Use in=construction_members.fa out=fresh_profiles.tsv with no sampling, overwrite or append");
		c.fasta=parser.in1;
		c.output=parser.out1;
		return c;
	}

	private static void validate(final Config c){
		require(c!=null && c.fasta!=null && c.references!=null && c.assignments!=null
			&& c.output!=null && c.models>0 && c.members>0,
			"Declare the complete frozen construction join, reference and positive expected model/member counts");
		if(c.selected){
			require(!Shared.AMINO_IN && !Read.IUPAC_TO_N && !Read.U_TO_T && !Read.DOT_DASH_X_TO_N && !Read.LOWER_CASE_TO_N,
				"Selected roster voting preserves actual nucleotide IUPAC content without lossy input conversion");
			require(c.inputMembers>=c.members && c.ledger==null && c.census==null && c.summary==null && c.paddingBlacklist!=null
				&& (c.partitionBins==512 || c.partitionBins==1024) && !c.dual,
				"selected=t uses the new roster FASTA/assignments only, explicit expectedinput and partitionbins512/1024, and explicit paddingblacklist=none or a scoped table");
			if("none".equals(c.paddingBlacklist)){require(c.paddingFit==null, "An empty padding set has no learned fit scope");}
			else{require(c.paddingFit!=null && !c.paddingFit.isEmpty() && Tools.testInputFiles(false, true, c.paddingBlacklist), "A learned padding set requires its actual table and explicit scope");}
		}else{
			require(c.ledger!=null && c.census!=null && c.summary!=null && c.paddingBlacklist==null && c.paddingFit==null && c.inputMembers<0,
				"Legacy layouts retain their original full construction join; roster-only flags require selected=t");
		}
		require((c.partitionBins==0 || c.partitionBins==512 || c.partitionBins==1024)
			&& !(c.dual && c.partitionBins!=0),
			"Use partitionbins=512or1024 for the selected compact design; dual=t is a separate historical geometry");
		require(Tools.testInputFiles(false, true, c.fasta, c.references, c.ledger, c.assignments, c.census, c.summary),
			"Every frozen construction/reference input must be readable");
		require(Tools.testOutputFiles(false, false, false, c.output), "Profile output must be fresh and writable");
		require(Tools.testForDuplicateFiles(true, c.fasta, c.references, c.ledger, c.assignments, c.census, c.summary, c.output),
			"Construction/reference inputs and profile output must have distinct paths");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && DynamicDemiLog.exponentBits()==4 && !DynamicDemiLog.blacklistExists(),
			"Profile construction requires assertions, fixed exponent4 and no DDL blacklist");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Processing         ----------------*/
	/*--------------------------------------------------------------*/

	public static void build(final Config c){
		validate(c);
		final boolean previousTrim=Shared.TRIM_READ_DESCRIPTION;
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		try{
			Shared.TRIM_READ_DESCRIPTION=false;
			NcrnaDdlSavedData.sequences(c.references, read->{
				names.add(read.id);
				sequences.add(read.bases);
			});
			require(names.size()==c.models, "The complete actual reference roster must match expectedmodels");
			final String[] ids=names.toArray(new String[names.size()]);
			final byte[][] bases=sequences.toArray(new byte[sequences.size()][]);
			final NcrnaDdlMembers members=c.selected
				? NcrnaDdlMembers.loadRoster(c.fasta, c.assignments, ids, bases, c.inputMembers, c.members)
				: NcrnaDdlMembers.load(c.fasta, c.ledger, c.assignments, c.census, c.summary, ids, bases);
			require(members.keptMembers()==c.members, "The actual original kept-voter population must match expectedmembers");
			if(c.selected){
				final LongHashSet blacklist="none".equals(c.paddingBlacklist) ? new LongHashSet(4)
					: NcrnaDdlSketchBuilder.readBlacklist(c.paddingBlacklist, c.paddingFit);
				final NcrnaDdlSketchBuilder builder=new NcrnaDdlSketchBuilder(c.partitionBins, 12345, blacklist);
				final DynamicDemiLog voter=builder.newSketch();
				final NcrnaDdlSketchBuilder.Prefix prefix=new NcrnaDdlSketchBuilder.Prefix(c.partitionBins);
				final NcrnaDdlLayoutVotes votes=new NcrnaDdlLayoutVotes(c.partitionBins);
				final char[][] keys=new char[c.models][];
				final int[] counts=new int[c.models];
				for(int model=0; model<c.models; model++){
					votes.reset();
					for(final byte[] member:members.members(model)){
						builder.fill(member, voter);
						votes.add(voter.maxArray());
					}
					builder.fillPrefix(bases[model], prefix);
					final NcrnaDdlLayoutVotes.Result result=votes.finish(prefix, 100, true);
					keys[model]=result.keys;
					counts[model]=result.members;
				}
				NcrnaDdlProfileResource.writeSelected(c.output, keys, counts, ids, bases, c.partitionBins, blacklist);
			}else if(c.partitionBins!=0){
				final NcrnaDdlProfile[] profiles=NcrnaDdlModelProfiles.buildPartition(members, bases, 7, c.partitionBins);
				NcrnaDdlProfileResource.writeCompact(c.output, profiles, ids, bases, c.partitionBins);
			}else{
				final NcrnaDdlProfile[] profiles=c.dual
					? NcrnaDdlModelProfiles.buildDual(members, 9, 10, NcrnaDdlProfileResource.DUAL_BUCKETS, true)
					: NcrnaDdlModelProfiles.build(members, NcrnaDdlProfileResource.K, NcrnaDdlProfileResource.BUCKETS, false, true);
				NcrnaDdlProfileResource.write(c.output, profiles, ids, bases, c.dual);
			}
		}finally{Shared.TRIM_READ_DESCRIPTION=previousTrim;}
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String fasta, references, ledger, assignments, census, summary, output, paddingBlacklist, paddingFit;
		public int models=508, members=85726, partitionBins, inputMembers=-1;
		public boolean dual, selected;
	}
}
