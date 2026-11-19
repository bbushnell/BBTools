package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import map.LongHashSet;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GeneCaller;
import prok.FamilyResourceResolver;
import prok.NcrnaBoundaryRefiner;
import prok.NcrnaBoundaryTrace;
import prok.NcrnaFamily;
import prok.NcrnaFamilyLoader;
import prok.NcrnaReportedOrf;
import prok.NcrnaScavenger;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.Orf;
import shared.Shared;

/** Literal post-trim selection/report separation, plus both native allocation routes.
 * Mock proposals are deliberately independent of any trained network. No biological
 * accuracy, stochastic fitting or throughput claim is made by this finite fixture.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDeferredReportTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, NcrnaDeferredReportTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String familyRoot=null;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("familyroot")){familyRoot=b;}
				else{check(parser.parse(arg, a, b), "Unknown deferred-report fixture argument: "+arg);}
			}
			check(familyRoot!=null && parser.out1!=null && parser.in1==null && parser.in2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1 && new File(parser.out1).mkdir(), "Supply out=fresh_fixture_directory");
			for(int strand=0; strand<2; strand++){
				for(int mode=0; mode<4; mode++){seam(new File(parser.out1), strand, mode);}
				allocation(false, strand);
				allocation(true, strand);
				fixed(strand, 1);
				fixed(strand, 3);
			}
			packageCases(new File(parser.out1), new File(familyRoot));
			pp.outstream.println("NCRNA_DEFERRED_REPORT_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** A physical start remains fixed; two physical ends bypass refinement altogether. */
	private static void fixed(final int strand, final int mask){
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		final NcrnaFamily family=family(false);
		family.maxLen=100;
		family.trimAlignmentExtent=false;
		family.boundaryOnRawEndpoints=true;
		final Proposal proposal=new Proposal(false);
		family.boundaryRefinerFactory=ignored->proposal;
		final TrimSpy caller=new TrimSpy(family, false);
		GeneCaller.applyNcrnaFamilyControls(caller, family);
		final int stop=mask==3 ? 239 : 80;
		final NcrnaReportedOrf candidate=new NcrnaReportedOrf("literal", 0, stop, strand, bases);
		candidate.ncrnaFamily=family.name;
		candidate.trnaModel=NAME;
		caller.finishAcceptedBoundary(candidate, bases, 0, 0, 239, 0, .8f, stop+1, .8f, true, null, null, mask, null);
		check(candidate.start==0 && candidate.stop==stop && proposal.calls==(mask==3 ? 0 : 1),
			"Physical masks preserve baseline selection and fully fixed calls never invoke");
		candidate.applyReportedBounds();
		check(candidate.start==0 && candidate.stop==(mask==3 ? 239 : 83), "Only the available report end may move");
	}

	/** Three tiny symlink packages reuse existing full-error payloads; no fitting or source copies. */
	private static void packageCases(final File root, final File original) throws Exception{
		for(int policy=0; policy<3; policy++){
			final File group=new File(root, "loaded_"+policy);
			check(group.mkdir(), "Each package-propagation case has a fresh data directory");
			ReadWrite.copyFile(new File(original, "manifest.tsv").getPath(), new File(group, "manifest.tsv").getPath());
			Files.createSymbolicLink(new File(group, "seeds.fa.gz").toPath(), new File(original, "seeds.fa.gz").toPath());
			for(final String member:new String[]{"RF00174", "RF01482", "RF01689"}){
				final File source=new File(original, member), dest=new File(group, member);
				if(!member.equals("RF00174")){
					Files.createSymbolicLink(dest.toPath(), source.toPath());
					continue;
				}
				check(dest.mkdir(), "The selected member needs a separate configuration");
				final File[] files=source.listFiles();
				check(files!=null && files.length>0, "The existing member package must expose its payload bindings");
				for(final File file:files){
					if(!file.getName().equals("positional.config")){Files.createSymbolicLink(new File(dest, file.getName()).toPath(), file.toPath());}
				}
				try(Input in=new Input(new File(source, "positional.config").getPath());
					Output out=new Output(new File(dest, "positional.config").getPath())){
					for(byte[] row; (row=in.reader.nextLine())!=null;){out.writer.println(row);}
					if(policy>0){out.writer.println("deferreport="+(policy==2));}
				}
			}
			final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
			options.resources.root=group.getPath();
			options.endpoints=false;
			final NcrnaFamily loaded=NcrnaFamilyLoader.load("cobalamin_rf00174", options, new FamilyResourceResolver());
			check(loaded.deferBoundaryReport==(policy==2) && loaded.boundaryRefinerFactory!=null,
				"Absent/false/true deferreport must survive real full-error family loading without changing factory selection");
			final NcrnaScavenger caller=new NcrnaScavenger(loaded.library, loaded.models, loaded.modelNames,
				loaded.kmerSet, loaded.kLong, loaded.minLen, loaded.windowPad, loaded.indexK, loaded.indexTopN,
				loaded.adaptive, loaded.adaptFloor, loaded.adaptTopFrac, loaded.adaptQFrac, loaded.fixedMinHits,
				loaded.scoreA, loaded.scoreB, loaded.idPass, loaded.idBorderline);
			GeneCaller.applyNcrnaFamilyControls(caller, loaded);
			check(caller.deferBoundaryReport==(policy==2), "The exact loaded policy must reach the actual scavenger controls");
		}
	}

	/** Modes: successful trim, no invocation, raw invocation, invalid report length. */
	private static void seam(final File root, final int strand, final int mode){
		final byte[] bases=new byte[240];
		Arrays.fill(bases, (byte)'A');
		final NcrnaFamily family=family(false);
		family.maxLen=100;
		family.boundaryOnRawEndpoints=mode==2;
		final Proposal proposal=new Proposal(mode==3);
		family.boundaryRefinerFactory=ignored->proposal;
		final TrimSpy caller=new TrimSpy(family, mode==0 || mode==3);
		GeneCaller.applyNcrnaFamilyControls(caller, family);
		// The worker already owns the literal spy. Do not advertise this mock as
		// a supported production positional resource to the strict trace loader.
		family.boundaryRefinerFactory=null;
		final String path=new File(root, "seam_"+strand+"_"+mode+".tsv").getPath();
		try(Output output=new Output(path)){
			output.writer.print(NcrnaBoundaryTrace.INVOCATION_HEADER);
			final NcrnaBoundaryTrace.Writer writer=new NcrnaBoundaryTrace.Writer(output.writer, family.name, true);
			writer.input(0, "literal.fa");
			final NcrnaBoundaryTrace.Worker worker=writer.worker(0, 0);
			final NcrnaBoundaryTrace.Family trace=worker.family(family);
			caller.setBoundaryTrace(trace);
			worker.begin("literal", bases.length, 0);
			try{
				final NcrnaReportedOrf candidate=new NcrnaReportedOrf("literal", 75, 145, strand, bases);
				candidate.ncrnaFamily=family.name;
				candidate.trnaModel=NAME;
				candidate.orfScore=42;
				final boolean trimmed=caller.finishAcceptedBoundary(candidate, bases, 0, 0, 239, 0, .8f, 71, .8f, true, null, null);
				final int start=trimmed ? 80 : 75, stop=trimmed ? 139 : 145;
				check(trimmed==(mode==0 || mode==3) && candidate.start==start && candidate.stop==stop
					&& candidate.orfScore==42, "Acceptance and path must retain the literal ordinary trim result and cached score");
				check(proposal.calls==(mode==1 ? 0 : 1), "Raw eligibility controls actual invocation without changing selection geometry");
				final int reportStart=mode==0 || mode==2 ? start-2 : start;
				final int reportStop=mode==0 || mode==2 ? stop+3 : stop;
				check(candidate.reportedSenseStart()==reportStart && candidate.reportedSenseStop()==reportStop,
					"No invocation saves the trimmed report; invalid length falls back to the same accepted baseline");
				trace.returned(candidate);
				if(strand==1){candidate.flip();}
				final ArrayList<Orf> selected=new ArrayList<Orf>();
				selected.add(candidate);
				NcrnaReportedOrf.applyTo(selected);
				check(candidate.start==(strand==0 ? reportStart : 239-reportStop)
					&& candidate.stop==(strand==0 ? reportStop : 239-reportStart), "Report application must use flip parity on both strands");
				check(!candidate.applyReportedBounds(), "Repeated final reporting is idempotent");
				worker.finish(selected);
			}finally{worker.clear();}
		}
		try(Input input=new Input(path)){
			int candidates=0, totals=0;
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length>1 && line[0]=='C'){candidates++;}
				if(line.length>1 && line[0]=='T'){totals++;}
			}
			check(candidates==1 && totals==1, "Each completed seam retains one candidate and one contig conservation record");
		}
	}

	/** Both primary/rescue allocation and reused-Quantum allocation must produce the subtype. */
	private static void allocation(final boolean reused, final int strand){
		final NcrnaFamily family=family(reused);
		family.trimAlignmentExtent=false;
		family.boundaryOnRawEndpoints=true;
		family.boundaryRefinerFactory=ignored->new Proposal(false);
		final NcrnaScavenger caller=new NcrnaScavenger(family.library, null, family.modelNames,
			family.kmerSet, 17, 16, 64, 7, 1, false, 1, 0, 0, 1, 0, 1, .6f, .6f);
		GeneCaller.applyNcrnaFamilyControls(caller, family);
		final byte[] bases=new byte[224];
		Arrays.fill(bases, (byte)'N');
		System.arraycopy(family.library[0], 0, bases, 80, 64);
		final ArrayList<Orf> calls=caller.scavenge("literal", bases, strand, new ArrayList<int[]>());
		check(calls.size()==1 && calls.get(0) instanceof NcrnaReportedOrf
			&& calls.get(0).start==80 && calls.get(0).stop==143,
			"Native allocation and claims retain the exact64-base accepted gene; reused="+reused+", strand="+strand);
		final NcrnaReportedOrf result=(NcrnaReportedOrf)calls.get(0);
		check(result.reportedSenseStart()==78 && result.reportedSenseStop()==146,
			"Both allocation routes retain the independent -2/+3 proposal without changing physical presence");
	}

	private static NcrnaFamily family(final boolean reused){
		final byte[] gene="GATCCGTACGAGTCCATGCTAGACCTGATCGAGCTTAGGCATACGTCATGGTACCGATCAGTGC".getBytes(StandardCharsets.US_ASCII);
		long key=0;
		for(int i=0; i<17; i++){
			final byte base=gene[i];
			key=(key<<2)|(base=='A' ? 0 : base=='C' ? 1 : base=='G' ? 2 : 3);
		}
		final LongHashSet seeds=new LongHashSet();
		seeds.add(key);
		final NcrnaFamily result=new NcrnaFamily("cobalamin_rf00174", new byte[][]{gene}, null,
			new String[]{NAME}, seeds, 17, 16, 64, 7, 1, false, 1, 0, 0, 1, 0, 1, .6f, .6f);
		result.maxLen=reused ? 100 : Integer.MAX_VALUE;
		result.scavengePass2=false;
		result.reuseConsensusAlignment=reused;
		result.deferBoundaryReport=true;
		return result;
	}

	private static final class Proposal implements NcrnaBoundaryRefiner{
		Proposal(final boolean invalid_){invalid=invalid_;}
		@Override
		public boolean applies(final int model){return model==0;}
		@Override
		public boolean afterAlignmentTrim(){return true;}
		@Override
		public boolean refine(final Orf orf, final byte[] bases, final int model, final float gc,
				final float identity, final boolean quantum, final int fixed, final int min, final int max){
			check(model==0 && (fixed&~3)==0, "The literal seam requires its actual selected model and physical mask");
			calls++;
			if((fixed&1)==0){orf.start=invalid ? 0 : orf.start-2;}
			if((fixed&2)==0){orf.stop=invalid ? bases.length-1 : orf.stop+3;}
			return true;
		}
		private final boolean invalid;
		int calls;
	}

	private static final class TrimSpy extends NcrnaScavenger{
		TrimSpy(final NcrnaFamily family, final boolean succeeds_){
			super(family.library, null, family.modelNames, family.kmerSet, 17, 16, 64);
			succeeds=succeeds_;
		}
		@Override
		public boolean trimToAlignmentExtent(final Orf candidate, final byte[] bases, final int model,
				final int start, final int stop, final float identity){
			check(model==0 && identity==.8f, "Trim spy requires the actual declared selected model and identity");
			if(succeeds){
				candidate.start=80;
				candidate.stop=139;
			}
			return succeeds;
		}
		private final boolean succeeds;
	}

	private static void check(final boolean condition, final String why){if(!condition){throw new AssertionError(why);}}
	private static final String NAME="cobalamin_rf00174_literal";
}
