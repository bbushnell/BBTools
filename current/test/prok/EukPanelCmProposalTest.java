package test.prok;

import prok.*;

import java.io.File;
import fileIO.ByteFile;
import parse.LineParser1;

/** Tests preservation, physical-vs-biological distinctions and fail-loud joins.
 * @author Ganyu */
public final class EukPanelCmProposalTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="),"Expected outdir=fresh_directory");
		final String root=args[0].substring(7);require(new File(root).mkdir(),"Fresh fixture directory required");
		final String fixture=root+"/evidence";
		EukPanelCmEvidenceTest.main(new String[]{"outdir="+fixture});
		run(fixture,fixture+"/pass.hit_overlaps.tsv.gz",fixture+"/pass.summary.tsv",root+"/proposal");
		verify(fixture,root+"/proposal");
		require("signal recognition particle RNA".equals(EukPanelCmProposal.gffAttribute("ID=x;product=signal recognition particle RNA;Dbxref=RFAM:RF00017","product")),"GFF product spaces must survive until the semicolon");
		require(EukPanelCmProposal.gffAttribute("Note=product=SRP RNA;other_product=SRP RNA","product")==null,"Attribute name must match at an exact semicolon boundary");
		require("SRP%20RNA".equals(EukPanelCmProposal.gffAttribute("product=SRP%20RNA","product")),"Original encoded value must remain unchanged");
		boolean duplicate=false;try{EukPanelCmProposal.gffAttribute("product=one;product=two","product");}catch(IllegalArgumentException e){require(e.getMessage().contains("Duplicate GFF attribute"),"Wrong duplicate-key failure");duplicate=true;}
		require(duplicate,"Conflicting repeated GFF fields must fail rather than silently choose");
		for(String attrs:new String[]{"product;ID=x","ID=x;product="}){boolean malformed=false;try{EukPanelCmProposal.gffAttribute(attrs,"product");}catch(IllegalArgumentException e){require(e.getMessage().contains("Malformed"),"Wrong malformed-key failure");malformed=true;}require(malformed,"Malformed requested attribute is not an absent attribute");}
		require(EukPanelCmProposal.srpSubtype("Dbxref=RFAM:RF00169;product=signal recognition particle RNA",true).equals("SRP_small"),"Canonical small SRP RFAM");
		require(EukPanelCmProposal.srpSubtype("Dbxref=GeneID:1,RFAM:RF01857;product=SRP RNA",true).equals("SRP_large"),"Canonical large SRP RFAM");
		require(EukPanelCmProposal.srpSubtype("Dbxref=RFAM:RF00169,RFAM:RF01854;product=SRP RNA",true).equals("AMBIGUOUS_SRP_SUBTYPE"),"Competing SRP subtypes remain unresolved");
		require(EukPanelCmProposal.srpSubtype("Dbxref=RFAM:RF00017;product=signal recognition particle RNA",true).equals("OUTSIDE_CURRENT_SUBTYPE_MAP"),"Outside-map eukaryotic SRP remains visible");
		require(EukPanelCmProposal.srpSubtype("product=SRP RNA",true).equals("UNRESOLVED_SRP_SUBTYPE"),"Missing RFAM does not become guessed subtype");
		require(EukPanelCmProposal.srpSubtype("ncRNA_class=SRP_RNA;Dbxref=RFAM:RF00017",true).equals("OUTSIDE_CURRENT_SUBTYPE_MAP"),"Source class can identify outside-scope SRP without choosing target subtype");
		final String sourceSrp=read(fixture+"/pass.sources.tsv.gz")
			.replace("ID=f2\t","ID=f2;Dbxref=RFAM:RF00017;product=signal recognition particle RNA\t")
			.replace("ID=f3\t","ID=f3;Dbxref=RFAM:RF00169,RFAM:RF01854;product=SRP RNA\t")
			.replace("ID=f4\t","ID=f4;product=SRP RNA\t");
		EukPanelCmEvidenceTest.write(root+"/sources_srp.tsv",sourceSrp);
		run(fixture,fixture+"/pass.hit_overlaps.tsv.gz",fixture+"/pass.summary.tsv",root+"/srp",root+"/sources_srp.tsv");
		verifySrpRows(root+"/srp",sourceSrp);
		final String pair=read(fixture+"/pass.hit_overlaps.tsv.gz");final int newline=pair.indexOf('\n');
		EukPanelCmEvidenceTest.write(root+"/duplicate_pairs.tsv",pair+pair.substring(newline+1));
		expectFailure(fixture,root+"/duplicate_pairs.tsv",fixture+"/pass.summary.tsv",root+"/bad_pair","Repeated/reversed pair");
		final String summary=read(fixture+"/pass.summary.tsv");final StringBuilder missing=new StringBuilder();
		for(String line:summary.split("\n")){if(!line.contains("\tRF00005\t")){missing.append(line).append('\n');}}
		EukPanelCmEvidenceTest.write(root+"/missing_zero.tsv",missing.toString());
		expectFailure(fixture,fixture+"/pass.hit_overlaps.tsv.gz",root+"/missing_zero.tsv",root+"/bad_zero","Missing zero-hit");
		System.out.println("PASS original hit/source bytes, both coordinate systems, cross-family flags, unresolved edge status, zero-hit denominators, no admission, SRP RFAM subtype and duplicate/missing-scan failures");
	}
	static void run(String root,String pairs,String summary,String out){
		run(root,pairs,summary,out,root+"/pass.sources.tsv.gz");
	}
	static void run(String root,String pairs,String summary,String out,String features){
		EukPanelCmProposal.main(new String[]{"inputs="+root+"/inputs.tsv","sources="+root+"/sources.tsv","contigs="+root+"/contigs.tsv","contigsources="+root+"/sources.tsv",
			"hits="+root+"/pass.hits.tsv.gz","pairs="+pairs,"summary="+summary,"features="+features,"associations="+root+"/pass.associations.tsv.gz",
			"expected=1","expectedmodels=3","expectedruns=1","expectedhits=5","expectedsources=6","expectedassociations=3","outprefix="+out});
	}
	static void verify(String fixture,String out){
		final ByteFile old=RrnaResourceIO.open(fixture+"/pass.hits.tsv.gz"),now=RrnaResourceIO.open(out+".overlay.tsv.gz");final LineParser1 p=new LineParser1('\t');int n=0;
		try{old.nextLine();p.set(now.nextLine());require(p.terms()==40,"Overlay has29original plus11descriptive fields");
			for(byte[] line;(line=now.nextLine())!=null;){final byte[] prior=old.nextLine();prefix(prior,line);p.set(line);
				require(p.terms()==40 && p.termEquals("REVIEW_ONLY",39) && p.termEquals("test",29),"No hit can be silently admitted");
				if(n<2){require(p.parseLong(38)==1 && p.parseLong(37)==0 && p.parseLong(36)==0,"Cross-family pair must remain distinct from same-family support");}
				else{require(p.parseLong(38)==0,"Unrelated or unresolved hit must not gain ambiguity edges");}
				if(n>=3){require(p.termEquals("UNASSESSABLE",31),"Invalid or unresolved interval has no inferred contig-edge contact");}
				n++;
			}require(n==5 && old.nextLine()==null,"Every raw hit retained exactly once");
		}finally{require(!old.close() && !now.close(),"Overlay fixture read failed");}
		final String[] sourceRows=read(fixture+"/pass.sources.tsv.gz").split("\n");
		final ByteFile links=RrnaResourceIO.open(out+".source_links.tsv.gz");n=0;
		try{p.set(links.nextLine());require(p.terms()==38,"Both-coordinate source-link schema");for(byte[] line;(line=links.nextLine())!=null;){
			p.set(line);require(p.terms()==38 && p.termEquals("PHYSICAL_LINK_REVIEW_ONLY",37),"Physical associations are never unique confirmed genes");
			final int sourceLine=(int)p.parseLong(4);prefix(sourceRows[sourceLine-1].getBytes(java.nio.charset.StandardCharsets.UTF_8),line);
			if(n==2){require(p.parseLong(31)==60 && p.parseLong(32)==40 && p.termEquals("-",33),"Minus CM coordinates remain original orientation beside source40..60");}
			n++;
		}require(n==3,"All three positive physical links retained");}finally{require(!links.close(),"Link read failed");}
		final ByteFile missing=RrnaResourceIO.open(out+".source_review.tsv.gz");n=0;
		try{p.set(missing.nextLine());require(p.terms()==28,"Source review schema");for(byte[] line;(line=missing.nextLine())!=null;){p.set(line);require(p.terms()==28 && p.termEquals("UNRESOLVED_SOURCE_REVIEW_ONLY",27),"Source anomalies stay unresolved");n++;}require(n==4,"Four no-link RNA/anomaly rows retained");}finally{require(!missing.close(),"Source review read failed");}
		final String[] denom=read(out+".scan_denominators.tsv").split("\n");require(denom.length==4 && denom[3].contains("RF00005"),"Zero-hit model remains in denominator");
	}
	static void prefix(byte[] expected,byte[] actual){require(expected!=null && actual.length>expected.length && actual[expected.length]=='\t',"Original row must be a complete prefix");for(int i=0;i<expected.length;i++){require(expected[i]==actual[i],"Original evidence field changed");}}
	static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final StringBuilder b=new StringBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(new String(line,java.nio.charset.StandardCharsets.UTF_8)).append('\n');}}finally{require(!in.close(),"Fixture read failed");}return b.toString();}
	static void verifySrpRows(String out,String sources){
		final String[] originals=sources.split("\n");final LineParser1 p=new LineParser1('\t');
		final ByteFile in=RrnaResourceIO.open(out+".source_links.tsv.gz");int n=0;
		try{in.nextLine();for(byte[] line;(line=in.nextLine())!=null;){p.set(line);prefix(originals[(int)p.parseLong(4)-1].getBytes(java.nio.charset.StandardCharsets.UTF_8),line);
			require(p.termEquals(n<2?"OUTSIDE_CURRENT_SUBTYPE_MAP":"AMBIGUOUS_SRP_SUBTYPE",26),"Outside/ambiguous SRP source row was dropped or assigned a winner");n++;
		}require(n==3,"All SRP physical links retained");}finally{require(!in.close(),"SRP links read failed");}
		final ByteFile missing=RrnaResourceIO.open(out+".source_review.tsv.gz");n=0;boolean found=false;
		try{missing.nextLine();for(byte[] line;(line=missing.nextLine())!=null;){p.set(line);prefix(originals[(int)p.parseLong(4)-1].getBytes(java.nio.charset.StandardCharsets.UTF_8),line);
			if(p.parseLong(4)==4){require(p.termEquals("UNRESOLVED_SRP_SUBTYPE",26),"Unassociated SRP lacking RFAM remains unresolved");found=true;}n++;
		}require(n==4 && found,"No-link source row conservation after SRP classification");}finally{require(!missing.close(),"SRP review read failed");}
	}
	static void expectFailure(String root,String pairs,String summary,String out,String reason){boolean failed=false;try{run(root,pairs,summary,out);}catch(IllegalArgumentException e){require(e.getMessage().contains(reason),"Wrong failure: "+e);failed=true;}require(failed,"Expected failure: "+reason);}
	static void require(boolean ok,String why){EukPanelCmProposal.require(ok,why);}
}
