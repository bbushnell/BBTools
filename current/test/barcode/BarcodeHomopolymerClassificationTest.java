package test.barcode;

import java.util.Locale;

import barcode.Barcode;
import barcode.BarcodeStats;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks homopolymer return codes and their weighted/unique statistics buckets.
 * The existing domain requires a defined first base, at least two matching bases,
 * no differing letter bases and at most one nonletter. Raw suffix text is distinct
 * from a plain barcode name with separately stored tile metadata.
 * @author Brian Bushnell, Sucrose
 */
public final class BarcodeHomopolymerClassificationTest{

	/** Requires mode=baseline|fixed; baseline qualifies the incorrect A fallback. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, BarcodeHomopolymerClassificationTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Classification checks require -ea.");
			Shared.setThreads(1);
			BarcodeStats.calcEdist=false;
			final int fallback="baseline".equals(mode) ? 0 : -1;
			int cases=0;
			final String[] rejected={"A", "C", "G", "T", "A+", "C+", "AA--AA", "CC++CC",
					"GG--GG", "TT++TT", "AAAA1101", "CCCC1101", "GGGG1101", "TTTT1101"};
			for(final String name : rejected){checkCase(name, fallback, false); cases++;}
			final String[] accepted={"AA", "CC", "GG", "TT", "AA+AA", "CC-CC", "GG+GG", "TT-TT", "CC1"};
			final int[] codes={0, 1, 2, 3, 0, 1, 2, 3, 1};
			for(int i=0; i<accepted.length; i++){checkCase(accepted[i], codes[i], false); cases++;}
			checkCase("AC", -1, false); cases++;
			checkCase("CC-GG", -1, false); cases++;
			checkCase("NN", -1, true); cases++;
			checkCase("CN", -1, true); cases++;
			checkCase("CC", 1101, "CC1101", 1, false); cases++;
			System.out.println("BARCODE_HOMOPOLYMER_CLASSIFICATION mode="+mode+" cases="+cases+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Raw names retain all characters; no CLI normalization or tile convention is inferred. */
	private static void checkCase(final String name, final int code, final boolean undefined){
		checkCase(name, 0, name, code, undefined);
	}

	/** Uses literal codes and weight seven; fresh instances avoid accumulated or expected-name shortcuts. */
	private static void checkCase(final String name, final int tile, final String key,
			final int code, final boolean undefined){
		assert(!name.isEmpty()) : "Barcode.checkHomopolymer reads the first character; fixture names must be nonempty.";
		final Barcode barcode=new Barcode(name, 7, 1, tile);
		final int actual=barcode.checkHomopolymer();
		check(actual==code, name+" tile="+tile+": expected code="+code+", observed="+actual);
		final BarcodeStats stats=new BarcodeStats(0, 1, null);
		if(tile==0){stats.increment(name, 7);}
		else{stats.increment(name, 7, tile);}
		stats.calcStats();
		check(stats.totalCodes==7 && stats.totalCodesU==1 && stats.expectedCodes==0 && stats.expectedCodesU==0,
				name+": totals must be weight7/unique1 with no expected observations; observed="+
				stats.totalCodes+"/"+stats.totalCodesU+"/"+stats.expectedCodes+"/"+stats.expectedCodesU);
		check(stats.nCodes==(undefined ? 7 : 0) && stats.nCodesU==(undefined ? 1 : 0),
				name+": undefined weight/unique mismatch="+stats.nCodes+"/"+stats.nCodesU);
		for(int base=0; base<4; base++){
			final boolean polymer=!undefined && base==code;
			check(stats.polymerArray[base]==(polymer ? 7 : 0) && stats.polymerArrayU[base]==(polymer ? 1 : 0),
					name+": base="+base+", expected weight/unique="+(polymer ? "7/1" : "0/0")+
					", observed="+stats.polymerArray[base]+"/"+stats.polymerArrayU[base]);
		}
		final Barcode retained=stats.codeMap.get(key);
		check(stats.codeMap.size()==1 && retained!=null && retained.name.equals(name) && retained.tile==tile && retained.count()==7,
				name+": original name, tile, key and observation weight must be retained; expected key="+key+", tile="+tile);
		System.out.println("name="+name+" tile="+tile+" code="+actual+" undefined="+undefined+" weight=7 unique=1");
	}

	/** Keeps literal oracle failures visible independently of the assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private BarcodeHomopolymerClassificationTest(){}
}
