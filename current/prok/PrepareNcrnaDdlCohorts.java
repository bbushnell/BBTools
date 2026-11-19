package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlIdentity.Locus;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;

/** Freezes identities from saved captures without constructing sketch profiles.
 * Uses A_SampleByteFile's PreParser/Parser and native byte-stream lifecycle.
 * Optional accepted prior cohorts supply identical-source values without alignment.
 * @author Brian, Raiden
 */
public final class PrepareNcrnaDdlCohorts{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, PrepareNcrnaDdlCohorts.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String panels=null, selection=null, references=null, previousPanels=null, previousSelection=null, previous=null;
			int models=508, groups=0, previousGroups=0;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("panels")){panels=b;}
				else if(a.equals("selection")){selection=b;}
				else if(a.equals("ref")){references=b;}
				else if(a.equals("previouspanels")){previousPanels=b;}
				else if(a.equals("previousselection")){previousSelection=b;}
				else if(a.equals("previous")){previous=b;}
				else if(a.equals("expectedmodels")){models=Parse.parseIntKMG(b);}
				else if(a.equals("expectedgroups")){groups=Parse.parseIntKMG(b);}
				else if(a.equals("previousgroups")){previousGroups=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown cohort-preparation parameter: "+arg);}
			}
			require(panels!=null && selection!=null && references!=null && models>0 && groups>=0
				&& parser.out1!=null && parser.out2==null && parser.in1==null && parser.in2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use complete saved panels/selection, ref and out=fresh_prefix; no sampling, overwrite or append");
			require(previous==null ? previousPanels==null && previousSelection==null && previousGroups==0
				: previousPanels!=null && previousSelection!=null && previousGroups>0,
				"Identity reuse requires previous=accepted_prefix, previouspanels, previousselection and previousgroups together");
			RrnaResourceIO.fresh(parser.out1, ".identity.tsv", ".identity_bins.tsv", ".groups.tsv", ".preparation.tsv");
			final boolean trim=Shared.TRIM_READ_DESCRIPTION;
			final int zip=ReadWrite.MAX_ZIP_THREADS();
			final boolean zipSet=ReadWrite.SET_ZIP_THREADS;
			try{
				Shared.TRIM_READ_DESCRIPTION=false;
				ReadWrite.setZipThreads(1);
				prepare(panels, selection, references, models, groups==0 ? count(selection) : groups,
					previousPanels, previousSelection, previous, previousGroups, parser.out1);
			}finally{
				Shared.TRIM_READ_DESCRIPTION=trim;
				ReadWrite.setZipThreads(zip);
				ReadWrite.SET_ZIP_THREADS=zipSet;
			}
			pp.outstream.println("NCRNA_DDL_COHORT_PREPARATION_PASS_NO_SKETCH_SCORING");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Preparation         ----------------*/
	/*--------------------------------------------------------------*/

	private static void prepare(final String panels, final String selection, final String ref, final int models,
			final int groups, final String previousPanels, final String previousSelection,
			final String previous, final int previousGroups, final String output){
		final NcrnaDdlSavedData data=NcrnaDdlSavedData.load(panels, selection, models, groups);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> references=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(ref, read->{
			names.add(read.id);
			references.add(read.bases);
		});
		require(names.size()==models, "Preparation must retain the complete original consensus ordinal/header roster");
		final String[] headers=names.toArray(new String[names.size()]);
		final byte[][] bases=references.toArray(new byte[references.size()][]), labels=new byte[models][];
		for(int model=0; model<models; model++){labels[model]=headers[model].getBytes(StandardCharsets.UTF_8);}
		data.readNativeCounts(headers);
		final Locus[] values;
		if(previous==null){values=NcrnaDdlIdentity.measure(data, bases);}
		else{
			final NcrnaDdlSavedData old=NcrnaDdlSavedData.load(previousPanels, previousSelection, models, previousGroups);
			old.readNativeCounts(headers);
			final NcrnaDdlCohorts accepted=NcrnaDdlCohorts.read(previous, old, labels, bases);
			values=NcrnaDdlIdentity.extend(data, bases, old, accepted.identity);
		}
		final NcrnaDdlCohorts cohorts=new NcrnaDdlCohorts(data, values);
		cohorts.write(output, data, labels);
		int reused=0, measured=0, unavailable=0;
		for(final Locus value:values){
			if(value.reused()){reused++;}
			else if(value.status().equals(NcrnaDdlIdentity.MEASURED)){measured++;}
			else{unavailable++;}
		}
		assert(reused+measured+unavailable==data.parents.size())
			: "Preparation accounts for each saved parent once, including reused missing values";
		try(Output out=new Output(output+".preparation.tsv")){
			out.writer.println("assembly_groups\ttarget_loci\treused_values\tnew_measured_identities\tnew_unavailable_values");
			out.writer.print(new ByteBuilder().append(groups).tab().append(values.length).tab().append(reused)
				.tab().append(measured).tab().append(unavailable).nl());
		}
	}

	/** A small complete inventory census; the loader validates every binding next. */
	private static int count(final String path){
		int count=0;
		try(Input in=new Input(path)){
			NcrnaDdlSavedData.header(in, NcrnaDdlSavedData.INVENTORY_HEADER);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){count++;}
		}
		require(count>0, "Cohort preparation cannot infer an empty assembly population");
		return count;
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
}
