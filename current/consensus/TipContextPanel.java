package consensus;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/**
 * Prepares a bounded experimental tip-polishing panel and audits its edits.
 * Truth coordinates select the diagnostic population only; no truth base is
 * supplied to context discovery, construction, merging, or the polishing guard.
 * @author Fischl
 */
public final class TipContextPanel {

	/** Uses native argument expansion and keeps diagnostics separate from artifacts. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, TipContextPanel.class, false);
		try{new TipContextPanel(pp.args).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Accepts explicit paths; no stage can overwrite an existing output artifact. */
	private TipContextPanel(final String[] args){
		assert(args!=null) : "PreParser provides the expanded panel arguments.";
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] pair=arg.split("=", 2);
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("mode")){mode=b;}
			else if(a.equals("in")){input=b;}
			else if(a.equals("targets")){targets=b;}
			else if(a.equals("dir")){directory=b;}
			else if(a.equals("out")){output=b;}
			else if(a.equals("report")){report=b;}
			else if(a.equals("before")){before=b;}
			else if(a.equals("after")){after=b;}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		check(!parser.append, "Panel artifacts must not mix executions through append mode.");
		if("prepare".equals(mode)){
			check(input!=null && targets!=null && directory!=null, "prepare requires in=assembly targets=BED dir=NEW");
			check(Tools.testInputFiles(false, true, input, targets), "Unreadable preparation inputs.");
		}else if("merge".equals(mode)){
			check(input!=null && directory!=null && output!=null && report!=null, "merge requires in=assembly dir=panel out=FASTA report=edits");
			check(Tools.testForDuplicateFiles(true, input, output, report), "Outputs must not alias the original assembly.");
			check(Tools.testOutputFiles(false, false, false, output, report), "Merge outputs must be new.");
		}else if("grade".equals(mode)){
			check(input!=null && before!=null && after!=null && output!=null, "grade requires in=edits before=used_snps after=used_snps out=summary");
			check(Tools.testForDuplicateFiles(true, input, before, after, output), "Grade output must not alias its evidence.");
			check(Tools.testInputFiles(false, true, input, before, after), "Unreadable grading inputs.");
			check(Tools.testOutputFiles(false, false, false, output), "Grade output must be new.");
		}else{throw new IllegalArgumentException("mode must be prepare, merge or grade.");}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Each mode has a separate output boundary, allowing retained intermediate evidence. */
	private void process(){
		assert(mode!=null) : "The constructor validates a concrete panel mode.";
		if(mode.equals("prepare")){prepare();}
		else if(mode.equals("merge")){merge();}
		else{grade();}
	}

	/** Adds only first-level unique donors as controls, never a recursive neighbor expansion. */
	private void prepare(){
		final LinkedHashMap<String, Read> assembly=readFasta(input);
		final LinkedHashSet<String> selected=new LinkedHashSet<String>();
		final LineParser1 parser=new LineParser1('\t');
		for(byte[] line : ByteFile.toLines(targets)){
			if(line.length==0 || line[0]=='#'){continue;}
			parser.set(line);
			final String name=parser.parseString(0);
			check(assembly.containsKey(name), "Target not present in original assembly: "+name);
			selected.add(name);
		}
		check(!selected.isEmpty(), "An empty selected population cannot validate a polisher.");
		final LinkedHashSet<String> all=new LinkedHashSet<String>(selected);
		final LinkedHashMap<String, Match[]> cache=new LinkedHashMap<String, Match[]>();
		for(String name : selected){
			final Match[] matches=matches(assembly.get(name), assembly);
			cache.put(name, matches);
			for(Match match : matches){if(match.set.count==1){all.add(match.donor);}}
		}
		check(new File(directory).mkdir(), "Preparation directory must be new: "+directory);
		final ByteBuilder table=new ByteBuilder("id\tname\trole\tlength\tleft_offset\tright_offset\n");
		final ByteBuilder candidates=new ByteBuilder("id\tside\tcandidates\toverlap\tdonor\n");
		int index=0, contexts=0;
		for(String name : all){
			final Read target=assembly.get(name);
			Match[] matches=cache.get(name);
			if(matches==null){matches=matches(target, assembly);}
			final String id=String.format(Locale.ROOT, "case%03d", ++index), dir=directory+"/"+id;
			check(new File(dir).mkdir(), "Case directory must be new: "+dir);
			writeFasta(dir+"/original.fa", target.id, target.bases);
			write(dir+"/target.bed", new ByteBuilder().append(target.id).tab().append(0).tab().append(target.length()).nl());
			final int[] offsets={-1, -1};
			for(int side=0; side<2; side++){
				final Match match=matches[side];
				candidates.append(id).tab().append(side==0 ? "left" : "right").tab().append(match.set.count).tab();
				if(match.set.count==1){
					final TipContextGuard.Candidate candidate=match.set.only;
					final byte[] context=context(target.bases, candidate, side==0);
					offsets[side]=side==0 ? context.length-target.length() : 0;
					writeFasta(dir+"/"+(side==0 ? "left" : "right")+".context.fa", "context", context);
					candidates.append(candidate.overlap).tab().append(match.donor).nl();
					contexts++;
				}else{candidates.append(-1).tab().append(".").nl();}
			}
			table.append(id).tab().append(name).tab().append(selected.contains(name) ? "selected" : "donor_control")
					.tab().append(target.length()).tab().append(offsets[0]).tab().append(offsets[1]).nl();
		}
		write(directory+"/targets.tsv", table);
		write(directory+"/candidates.tsv", candidates);
		System.out.println("TIP_PANEL_PREPARED selected="+selected.size()+" cases="+all.size()+" contexts="+contexts);
	}

	/** Combines only changes agreeing with at least one end and opposed by no other proposed change. */
	private void merge(){
		final LinkedHashMap<String, Read> assembly=readFasta(input);
		final HashSet<String> seen=new HashSet<String>();
		final ByteBuilder edits=new ByteBuilder("contig\tposition1\toriginal\tguarded\trole\n");
		final LineParser1 parser=new LineParser1('\t');
		int changes=0, conflicts=0;
		for(byte[] line : ByteFile.toLines(directory+"/targets.tsv")){
			parser.set(line);
			if(parser.termEquals("id", 0)){continue;}
			check(parser.terms()==6, "Unexpected case manifest schema.");
			final String id=parser.parseString(0), name=parser.parseString(1), role=parser.parseString(2);
			check(seen.add(name) && assembly.containsKey(name), "Repeated or missing case target: "+name);
			final Read target=assembly.get(name), original=single(directory+"/"+id+"/original.fa");
			check(target.id.equals(original.id) && Arrays.equals(target.bases, original.bases), "Original assembly changed before merge: "+name);
			final byte[][] proposed=new byte[2][];
			for(int side=0; side<2; side++){
				if(parser.parseInt(4+side)<0){proposed[side]=original.bases; continue;}
				final Read guarded=single(directory+"/"+id+"/"+(side==0 ? "left" : "right")+".guarded.fa");
				check(guarded.id.equals(name) && guarded.length()==original.length(), "Guard changed target identity or length: "+name);
				proposed[side]=guarded.bases;
			}
			final byte[] result=combine(original.bases, proposed[0], proposed[1]);
			for(int pos=0; pos<result.length; pos++){
				if(proposed[0][pos]!=original.bases[pos] && proposed[1][pos]!=original.bases[pos]
						&& proposed[0][pos]!=proposed[1][pos]){conflicts++;}
				if(result[pos]==original.bases[pos]){continue;}
				edits.append(name).tab().append(pos+1).tab().append((char)original.bases[pos]).tab()
						.append((char)result[pos]).tab().append(role).nl();
				changes++;
			}
			target.bases=result;
		}
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, false);
		writer.start();
		final ByteBuilder record=new ByteBuilder();
		for(Read read : assembly.values()){
			record.clear();
			record.append('>').append(read.id).nl().append(read.bases).nl();
			writer.print(record);
		}
		check(!writer.poisonAndWait(), "Failed to finish the complete reconstructed assembly.");
		write(report, edits);
		System.out.println("TIP_PANEL_MERGED cases="+seen.size()+" substitutions="+changes+" conflicts="+conflicts);
	}

	/** Checks every changed QUAST mismatch locus against an actually emitted substitution. */
	private void grade(){
		final HashSet<String> old=snps(before), current=snps(after), edits=new HashSet<String>();
		final LineParser1 parser=new LineParser1('\t');
		int repaired=0, damaged=0, stillWrong=0, unclassified=0;
		for(byte[] line : ByteFile.toLines(input)){
			parser.set(line);
			if(parser.termEquals("contig", 0)){continue;}
			check(parser.terms()==5, "Unexpected edit ledger schema.");
			final String name=parser.parseString(0).replace(',', '_').replace('=', '_');
			final String key=name+"\t"+parser.parseInt(1);
			check(edits.add(key), "An assembly base must be edited at most once: "+key);
			if(old.contains(key)){
				if(current.contains(key)){stillWrong++;}else{repaired++;}
			}else if(current.contains(key)){damaged++;}
			else{unclassified++;}
		}
		for(String key : old){check(current.contains(key) || edits.contains(key), "An unedited old mismatch disappeared: "+key);}
		for(String key : current){check(old.contains(key) || edits.contains(key), "An unedited new mismatch appeared: "+key);}
		write(output, new ByteBuilder("before_errors\tafter_errors\tedits\trepaired\tnew_wrong\tchanged_still_wrong\tunclassified\n")
				.append(old.size()).tab().append(current.size()).tab().append(edits.size()).tab().append(repaired)
				.tab().append(damaged).tab().append(stillWrong).tab().append(unclassified).nl());
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Searches the frozen assembly, excluding only the target's own physical record. */
	private static Match[] matches(final Read target, final LinkedHashMap<String, Read> assembly){
		assert(target!=null && assembly.containsKey(target.id)) : "Context discovery requires a target from the candidate universe.";
		final Match[] result={new Match(), new Match()};
		for(int side=0; side<2; side++){
			final byte[] query=side==0 ? TipContextGuard.reverse(target.bases) : target.bases;
			final Match match=result[side];
			for(Read donor : assembly.values()){
				if(donor==target){continue;}
				final int prior=match.set.count;
				TipContextGuard.scan(query, donor.bases, 60, 200, match.set);
				TipContextGuard.scan(query, TipContextGuard.reverse(donor.bases), 60, 200, match.set);
				if(prior==0 && match.set.count==1){match.donor=donor.id;}
			}
		}
		return result;
	}

	/** Pads only one end, preserving every original target base in its original orientation. */
	static byte[] context(final byte[] target, final TipContextGuard.Candidate donor, final boolean left){
		assert(target!=null && donor!=null && donor.overlap<donor.bases.length) : "Context needs an extending donor, not a contained match.";
		final byte[] oriented=left ? TipContextGuard.reverse(target) : target;
		final byte[] joined=Arrays.copyOf(oriented, Math.addExact(target.length, donor.bases.length-donor.overlap));
		System.arraycopy(donor.bases, donor.overlap, joined, target.length, donor.bases.length-donor.overlap);
		return left ? TipContextGuard.reverse(joined) : joined;
	}

	/** Disagreeing changes abstain; unchanged sides do not veto an independent accepted edit. */
	static byte[] combine(final byte[] original, final byte[] left, final byte[] right){
		assert(original.length==left.length && original.length==right.length) : "Substitution-only merging requires coordinate parity.";
		final byte[] result=original.clone();
		for(int pos=0; pos<result.length; pos++){
			final byte a=left[pos], b=right[pos], old=original[pos];
			if(a==b){result[pos]=a;}
			else if(a==old){result[pos]=b;}
			else if(b==old){result[pos]=a;}
		}
		return result;
	}

	/** Uses native sequence streaming and retains whole named contigs, not hand-parsed FASTA. */
	private static LinkedHashMap<String, Read> readFasta(final String path){
		check(path!=null && new File(path).isFile(), "Missing sequence evidence: "+path);
		final FileFormat ff=FileFormat.testInput(path, FileFormat.FASTA, null, true, true);
		check(ff.fasta(), "Panel sequence inputs must be FASTA: "+path);
		final Streamer streamer=StreamerFactory.makeStreamer(ff, null, true, -1, false, true, 1);
		final LinkedHashMap<String, Read> result=new LinkedHashMap<String, Read>();
		streamer.start();
		try{
			for(ListNum<Read> batch=streamer.nextList(); batch!=null; batch=streamer.nextList()){
				for(Read read : batch){
					check(read.mate==null && read.bases!=null && result.put(read.id, read)==null, "Invalid/duplicate FASTA name: "+read.id);
				}
			}
		}finally{streamer.close();}
		check(!streamer.errorState() && !result.isEmpty(), "Empty or failed FASTA read: "+path);
		return result;
	}

	/** The case workflow never permits silent selection of a different record. */
	private static Read single(final String path){
		final LinkedHashMap<String, Read> reads=readFasta(path);
		check(reads.size()==1, "Expected exactly one target: "+path);
		return reads.values().iterator().next();
	}

	/** Extracts native QUAST substitution coordinates, keeping every distinct contig/base. */
	private static HashSet<String> snps(final String path){
		final HashSet<String> result=new HashSet<String>();
		final LineParser1 parser=new LineParser1('\t');
		for(byte[] line : ByteFile.toLines(path)){
			if(line.length==0){continue;}
			parser.set(line);
			check(parser.terms()>=6, "Malformed native QUAST used_snps row.");
			check(!parser.termEquals(".", 3) && !parser.termEquals(".", 4), "Indels require a different coordinate audit.");
			check(result.add(parser.parseString(1)+"\t"+parser.parseInt(5)), "Repeated QUAST contig/position.");
		}
		return result;
	}

	/** Emits one unwrapped record using the native buffered writer. */
	private static void writeFasta(final String path, final String name, final byte[] bases){
		assert(bases!=null && name!=null) : "Every panel FASTA must retain sequence and identity.";
		write(path, new ByteBuilder().append('>').append(name).nl().append(bases).nl());
	}

	/** Never overwrites or silently drops an asynchronous output error. */
	private static void write(final String path, final ByteBuilder contents){
		check(!new File(path).exists(), "Output already exists: "+path);
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		writer.print(contents);
		check(!writer.poisonAndWait(), "Could not finish panel artifact: "+path);
	}

	/** Evidence failures remain fatal even if an external caller disables assertions. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Tracks physical multiplicity separately from the first donor's identity. */
	private static final class Match {
		final TipContextGuard.CandidateSet set=new TipContextGuard.CandidateSet();
		String donor;
	}

	private String mode, input, targets, directory, output, report, before, after;
}
