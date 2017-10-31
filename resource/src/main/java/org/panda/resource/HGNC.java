package org.panda.resource;

import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.UndirectedGraph;

import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

/**
 * This class provides a mapping between HGNC IDs and approved gene symbols.
 *
 * @author Ozgun Babur
 */
public class HGNC extends FileServer
{
	private static HGNC instance;

	private Map<String, String> sym2id;
	private Map<String, String> sym2chr;
	private Map<String, String> id2sym;
	private Map<String, String> old2new;
	private Map<String, String> uniprot2sym;
	private Map<String, Set<String>> families = new HashMap<>();

	public static synchronized HGNC get()
	{
		if (instance == null) instance = new HGNC();
		return instance;
	}

	/**
	 * Gets the latest approved official symbol related to the given ID or symbol. If the parameter
	 * is ID, then it should start with "HGNC:".
	 * @param symbolOrID HGNC ID, symbol, or a previous symbol
	 * @return latest symbol
	 */
	public String getSymbol(String symbolOrID)
	{
		if (symbolOrID == null) return null;
		if (id2sym.containsKey(symbolOrID)) return id2sym.get(symbolOrID);
		else if (sym2id.containsKey(symbolOrID)) return symbolOrID;
		symbolOrID = symbolOrID.toUpperCase();
		if (old2new.containsKey(symbolOrID)) return old2new.get(symbolOrID);
		if (uniprot2sym.containsKey(symbolOrID)) return uniprot2sym.get(symbolOrID);
		return null;
	}

	public Set<String> getFamily(String name)
	{
		if (!families.containsKey(name)) return Collections.emptySet();
		return families.get(name);
	}

	public Set<String> getAllSymbols()
	{
		return sym2id.keySet();
	}

	public Set<String> getSymbolsOfChromosome(String chr)
	{
		Set<String> set = new HashSet<String>();
		for (String sym : getAllSymbols())
		{
			if (sym2chr.containsKey(sym))
			{
				String c = sym2chr.get(sym);
				String no = c.split("p")[0].split("q")[0];

				String desNo = chr.split("p")[0].split("q")[0];

				if (no.equals(desNo) && c.contains(chr))
				{
					set.add(sym);
				}
			}
		}
		return set;
	}

	public String getChromosomeLoc(String symbol)
	{
		return sym2chr.get(symbol);
	}

	/**
	 * If the given collection
	 */
	public void replaceRecognizedOlds(Collection<String> col, boolean removeUnrecognized)
	{
		Set<String> news = new HashSet<>();
		Set<String> discard = new HashSet<>();
		for (String s : col)
		{
			String n = getSymbol(s);
			if (n == null)
			{
				if (removeUnrecognized) discard.add(s);
			}
			else if (!n.equals(s))
			{
				news.add(n);
				discard.add(s);
			}
		}
		col.removeAll(discard);
		col.addAll(news);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"hgnc.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym" +
			"&col=gd_aliases&col=gd_pub_chrom_map&col=family.name&col=md_prot_id&status=Approved&status_opt=2&where=&" +
			"order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"};
	}

	@Override
	public boolean load() throws IOException
	{
		sym2id = new HashMap<>();
		sym2chr = new HashMap<>();
		id2sym = new HashMap<>();
		old2new = new HashMap<>();
		uniprot2sym = new HashMap<>();
		families = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(1).forEach(line -> {
			String[] token = line.split("\t");
			String sym = token[1];
			String id = token[0];
			sym2id.put(sym, id);
			id2sym.put(id, sym);

			if (token.length > 2)
			{
				for (String old : token[2].split(","))
				{
					old = old.trim().toUpperCase();
					if (old.isEmpty()) continue;
					old2new.put(old, sym);
				}
			}
			if (token.length > 3)
			{
				for (String synonym : token[3].split(","))
				{
					synonym = synonym.trim().toUpperCase();
					if (synonym.isEmpty()) continue;
					old2new.put(synonym, sym);
				}
			}
			if (token.length > 4)
			{
				sym2chr.put(sym, token[4]);
			}
			if (token.length > 5)
			{
				if (!families.containsKey(token[5]))
					families.put(token[5], new HashSet<>());

				families.get(token[5]).add(sym);
			}
			if (token.length > 6 && !token[6].isEmpty())
			{
				uniprot2sym.put(token[6], sym);
			}

		});
		return true;
	}

	public Graph getCompleteClique(boolean directed)
	{
		Graph graph = directed ? new DirectedGraph("HGNC complete clique", "clique-edge") :
			new UndirectedGraph("HGNC complete clique", "clique-edge");

		for (String s1 : sym2id.keySet())
		{
			for (String s2 : sym2id.keySet())
			{
				if (s1.equals(s2)) continue;
				if (!directed && s2.compareTo(s1) < 0) continue;
				graph.putRelation(s1, s2);
				if (directed) graph.putRelation(s2, s1);
			}
		}
		return graph;
	}

	public static void main(String[] args)
	{
//		HGNC hgnc = new HGNC();
//		System.out.println(hgnc.getSymbol("CDC42"));
//		System.out.println(hgnc.getChromosomeLoc("ATAD2"));
//		Set<String> set = hgnc.getSymbolsOfChromosome("8q24");
//		for (String sym : set)
//		{
//			System.out.println(sym + "\t" + hgnc.getChromosomeLoc(sym));
//		}
//
//		System.out.println(hgnc.getChromosomeLoc("MYC"));

		convertSymbols();
	}

	private static void convertSymbols()
	{
		List<String> list = Arrays.asList(("TAL1\n" +
			"AGBL4\n" +
			"ELAVL4\n" +
			"NEGR1\n" +
			"FPGT-TNNI3K\n" +
			"FUBP1\n" +
			"PTBP2\n" +
			"GNAT2\n" +
			"SEC16B\n" +
			"NAV1\n" +
			"TMEM18\n" +
			"ADCY3\n" +
			"KCNK3\n" +
			"LINC01122\n" +
			"EHBP1\n" +
			"LRP1B\n" +
			"FIGN\n" +
			"COBLL1-GRB14\n" +
			"UBE2E3\n" +
			"CREB1\n" +
			"ERBB4\n" +
			"USP37\n" +
			"LOC646736\n" +
			"IRS1a\n" +
			"RARB\n" +
			"FHIT\n" +
			"GBE1\n" +
			"CADM2\n" +
			"RASA2\n" +
			"ETV5\n" +
			"GNPDA2\n" +
			"SCARB2\n" +
			"SLC39A8\n" +
			"HHIP\n" +
			"POC5\n" +
			"GALNT10\n" +
			"C6orf106\n" +
			"TDRG1\n" +
			"TFAP2B\n" +
			"FOXO3\n" +
			"LOC285762\n" +
			"IFNGR1\n" +
			"PARK2\n" +
			"HIP1\n" +
			"PMS2L11\n" +
			"CALCR\n" +
			"ASB4\n" +
			"HNF4G\n" +
			"ZBTB10\n" +
			"RALYL\n" +
			"C9orf93\n" +
			"LINGO2\n" +
			"EPB41L4B\n" +
			"TLR4\n" +
			"LMX1B\n" +
			"GRID1\n" +
			"HIF1AN\n" +
			"NT5C2\n" +
			"TCF7L2\n" +
			"TRIM66\n" +
			"BDNF\n" +
			"HSD17B12\n" +
			"MTCH2\n" +
			"CADM1\n" +
			"BCDIN3D\n" +
			"CLIP1\n" +
			"MTIF3\n" +
			"OLFM4\n" +
			"MIR548X2\n" +
			"MIR548A2\n" +
			"SPRY2a\n" +
			"STXBP6\n" +
			"PRKD1\n" +
			"PRKD1\n" +
			"NRXN3\n" +
			"DMXL2\n" +
			"MAP2K5\n" +
			"LOC100287559\n" +
			"NLRC3\n" +
			"GPRC5B\n" +
			"SBK1\n" +
			"ATP2A1\n" +
			"INO80E\n" +
			"KAT8\n" +
			"CBLN1\n" +
			"FTO\n" +
			"SMG6\n" +
			"RABEP1\n" +
			"IGF2BP1\n" +
			"RPTOR\n" +
			"C18orf8\n" +
			"LOC284260\n" +
			"GRP\n" +
			"MC4R\n" +
			"PGPEP1\n" +
			"CRTC1\n" +
			"KCTD15\n" +
			"TOMM40\n" +
			"QPCTL\n" +
			"ZC3H4\n" +
			"ZFP64\n" +
			"ETS2\n" +
			"PICK1-PLA2G6\n").split("\n"));

		for (String text : list)
		{
			String sym = get().getSymbol(text);
			System.out.println(text + "\t" + sym);
		}
	}
}
