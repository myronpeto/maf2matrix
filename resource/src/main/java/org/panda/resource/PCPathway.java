package org.panda.resource;

import org.biopax.paxtools.controller.Cloner;
import org.biopax.paxtools.controller.Completer;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.*;
import org.biopax.paxtools.pattern.util.Blacklist;
import org.panda.utility.Kronometre;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.FishersExactTest;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Provides pathways and their contents in Pathway Commons. Also runs enrichment analyses.
 *
 * @author Ozgun Babur
 */
public class PCPathway extends FileServer
{
	private static PCPathway instance;
	private static final String FILE = "pcpathway.txt";

	private Map<String, Set<String>> gene2pathway;
	private Map<String, Set<String>> pathway2gene;
	private Map<String, Set<String>> chem2pathway;
	private Map<String, Set<String>> pathway2chem;
	private Map<String, String> pathway2name;
	private Map<String, String> pathway2resource;

	public static synchronized PCPathway get()
	{
		if (instance == null) instance = new PCPathway();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{FILE};
	}

	@Override
	public boolean load() throws IOException
	{
		pathway2gene = new HashMap<>();
		gene2pathway = new HashMap<>();
		chem2pathway = new HashMap<>();
		pathway2chem = new HashMap<>();
		pathway2name = new HashMap<>();
		pathway2resource = new HashMap<>();

		Set<Set<String>> groups = new HashSet<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(line -> line.split("\t")).forEach(token ->
		{
			pathway2name.put(token[0], token[1]);
			pathway2resource.put(token[0], token[2]);

			if (token.length > 3)
			{
				Set<String> group = new HashSet<>(Arrays.asList(token).subList(3, token.length));
				if (groups.contains(group)) return;
				groups.add(group);

				Set<String> chems = Arrays.asList(token).subList(3, token.length).stream()
					.filter(s -> s.startsWith("CHEBI:")).collect(Collectors.toSet());

				Set<String> genes = Arrays.asList(token).subList(3, token.length).stream()
					.filter(s -> !s.startsWith("CHEBI:")).collect(Collectors.toSet());

				pathway2gene.put(token[0], genes);
				pathway2chem.put(token[0], chems);

				genes.stream().filter(gene -> !gene2pathway.containsKey(gene))
					.forEach(gene -> gene2pathway.put(gene, new HashSet<>()));
				genes.stream().forEach(gene -> gene2pathway.get(gene).add(token[0]));

				chems.stream().filter(chem -> !chem2pathway.containsKey(chem))
					.forEach(chem -> chem2pathway.put(chem, new HashSet<>()));
				chems.stream().forEach(chem -> chem2pathway.get(chem).add(token[0]));
			}
		});

		return true;
	}

	public void addCustomPathway(String id, String name, Set<String> genes, Set<String> chems)
	{
		if (pathway2name.containsKey(id)) throw new IllegalArgumentException("Pathway already exists. ID = " + id);

		pathway2name.put(id, name);
		pathway2gene.put(id, genes);
		pathway2chem.put(id, chems);
		for (String gene : genes)
		{
			if (!gene2pathway.containsKey(gene)) gene2pathway.put(gene, new HashSet<>());
			gene2pathway.get(gene).add(id);
		}
		for (String chem : chems)
		{
			if (!chem2pathway.containsKey(chem)) chem2pathway.put(chem, new HashSet<>());
			chem2pathway.get(chem).add(id);
		}
		pathway2resource.put(id, "Custom");
	}

	public Set<String> getPathways(String gene)
	{
		if (gene2pathway.containsKey(gene)) return gene2pathway.get(gene);
		else return Collections.emptySet();
	}

	public Set<String> getGenes(String pathwayID)
	{
		if (pathway2gene.containsKey(pathwayID)) return pathway2gene.get(pathwayID);
		else return Collections.emptySet();
	}

	public Set<String> getChems(String pathwayID)
	{
		if (pathway2chem.containsKey(pathwayID)) return pathway2chem.get(pathwayID);
		else return Collections.emptySet();
	}

	public String getName(String id)
	{
		return pathway2name.get(id);
	}

	public String getResource(String id)
	{
		return pathway2resource.get(id);
	}

	public String getCoverageStr(String id, Set<String> genes)
	{
		if (!pathway2gene.containsKey(id)) return null;
		Set<String> set = new HashSet<String>(genes);
		set.retainAll(pathway2gene.get(id));
		return set.size() + "/" + pathway2gene.get(id).size();
	}

	/**
	 * Gets the enrichment pvals and pval limits.
	 * @param molecules query
	 * @param backgroundMols if there is any
	 * @return two maps, first is for pvals, second is for limits
	 */
	public Map<String, Double>[] getEnrichmentPvals(Collection<String> molecules,
		Collection<String> backgroundMols, int minMemberSize, int maxMemberSize)
	{
		if (backgroundMols != null && !backgroundMols.containsAll(molecules)) throw new IllegalArgumentException(
			"Background molecules have to contain all the selected molecules.");

		Set<String> mols = new HashSet<>(molecules);
		Set<String> genes = mols.stream().filter(entity -> !entity.startsWith("CHEBI:")).collect(Collectors.toSet());
		Set<String> chems = mols.stream().filter(entity -> entity.startsWith("CHEBI:")).collect(Collectors.toSet());

		Set<String> background = new HashSet<>();
		if (!genes.isEmpty()) background.addAll(gene2pathway.keySet());
		if (!chems.isEmpty()) background.addAll(chem2pathway.keySet());

		if (backgroundMols != null) background.retainAll(backgroundMols);

		if (!background.containsAll(mols))
		{
			Set<String> excess = new HashSet<>(mols);
			excess.removeAll(background);

			mols.removeAll(excess);
			genes.removeAll(excess);
			chems.removeAll(excess);
			System.out.println("Removed " + excess.size() + " unknown genes: " + excess);
			System.out.println("Using " + mols.size() + ": " + mols);
		}

		Map<String, Integer> selectionGeneCnt = count(genes, pathway2gene);
		Map<String, Integer> backgroundGeneCnt = count(background, pathway2gene);
		Map<String, Integer> selectionChemCnt = count(chems, pathway2chem);
		Map<String, Integer> backgroundChemCnt = count(background, pathway2chem);

		Map<String, Double> mapP = new HashMap<>();
		Map<String, Double> mapL = new HashMap<>();

		Set<Set<String>> memory = new HashSet<>();

		Stream.concat(selectionGeneCnt.keySet().stream(), selectionChemCnt.keySet().stream()).distinct().forEach(pathway ->
		{
			int pathwaySize = 0;
			if (!genes.isEmpty())
			{
				pathwaySize += pathway2gene.get(pathway).size();
				if (chems.isEmpty())
				{
					Set<String> cs = new HashSet<>(pathway2gene.get(pathway));
					if (memory.contains(cs)) return;
					else memory.add(cs);
				}
			}
			if (!chems.isEmpty())
			{
				pathwaySize += pathway2chem.get(pathway).size();
				if (genes.isEmpty())
				{
					Set<String> cs = new HashSet<>(pathway2chem.get(pathway));
					if (memory.contains(cs)) return;
					else memory.add(cs);
				}
			}

			if (pathwaySize < minMemberSize || pathwaySize > maxMemberSize) return;

			int size = background.size();
			int featuredOverall = 0;
			if (!genes.isEmpty()) featuredOverall += backgroundGeneCnt.get(pathway);
			if (!chems.isEmpty()) featuredOverall += backgroundChemCnt.get(pathway);
			int selected = mols.size();
			int featuredSelected = 0;
			if (!genes.isEmpty()) featuredSelected += selectionGeneCnt.get(pathway);
			if (!chems.isEmpty()) featuredSelected += selectionChemCnt.get(pathway);

			double pval = FishersExactTest.calcEnrichmentPval(size, featuredOverall, selected, featuredSelected);

			int maxPossibleHit = 0;
			if (!genes.isEmpty()) maxPossibleHit += Math.min(backgroundGeneCnt.get(pathway), genes.size());
			if (!chems.isEmpty()) maxPossibleHit += Math.min(backgroundChemCnt.get(pathway), chems.size());

			double limit = FishersExactTest.calcEnrichmentPval(size, featuredOverall, selected, maxPossibleHit);

			mapP.put(pathway, pval);
			mapL.put(pathway, limit);
		});

		return new Map[]{mapP, mapL};
	}

	private Map<String, Integer> count(Collection<String> mols, Map<String, Set<String>> pathway2X)
	{
		Map<String, Integer> cnt = new HashMap<>();

		for (String pathway : pathway2X.keySet())
		{
			Set<String> mems = new HashSet<>(pathway2X.get(pathway));
			mems.retainAll(mols);
			cnt.put(pathway, mems.size());
		}
		return cnt;
	}

	public List<String> getEnrichedPathways(Collection<String> molecules,
		Collection<String> background, double fdrThr)
	{
		return getEnrichedPathways(molecules, background, fdrThr, 3, 500);
	}

	public List<String> getEnrichedPathways(Collection<String> molecules,
		Collection<String> background, double fdrThr, int minMemberSize, int maxMemberSize)
	{
		Map<String, Double>[] map = getEnrichmentPvals(molecules, background, minMemberSize, maxMemberSize);
		if (fdrThr < 0)
		{
			fdrThr = FDR.decideBestFDR_BH(map[0], map[1]);
			System.out.println("fdrThr = " + fdrThr);
		}
		return FDR.select(map[0], map[1], fdrThr).stream().collect(Collectors.toList());
	}

	public Map<String, Double> getEnrichmentQvals(Collection<String> molecules,
		Collection<String> background, int minMemberSize, int maxMemberSize)
	{
		Map<String, Double>[] map = getEnrichmentPvals(molecules, background, minMemberSize, maxMemberSize);
		return FDR.getQVals(map[0], map[1]);
	}

	public void writeEnrichmentResults(Set<String> molecules, int minMemberSize, int maxMemberSize, String filename)
		throws IOException
	{
		writeEnrichmentResults(molecules, null, minMemberSize, maxMemberSize,filename);
	}

	public void writeEnrichmentResults(Set<String> molecules, Set<String> background, int minMemberSize, int maxMemberSize, String filename)
		throws IOException
	{
		final Map<String, Double>[] pvals = getEnrichmentPvals(molecules, background, minMemberSize, maxMemberSize);
		Map<String, Double> qvals = getEnrichmentQvals(molecules, background, minMemberSize, maxMemberSize);

		boolean hasGene = molecules.stream().anyMatch(s -> !s.startsWith("CHEBI:"));
		boolean hasChem = molecules.stream().anyMatch(s -> s.startsWith("CHEBI:"));

		List<String> ids = new ArrayList<>(qvals.keySet());
		Collections.sort(ids, (o1, o2) -> pvals[0].get(o1).compareTo(pvals[0].get(o2)));

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		writer.write("# Pathway name: Name of the pathway as given by the original pathway database.\n");
		writer.write("# Resource: The original database that this pathway was defined.\n");
		writer.write("# Pathway Commons ID: ID of the pathway in Pathway Commons database.\n");
		writer.write("# P-value: Enrichment p-value of the pathway calculated by Fisher's exact test.\n");
		writer.write("# Q-value: Estimated FDR (false discovery rate) if this p-value is used as cutoff threshold.\n");
		writer.write("# Hit size: Number of query genes that overlaps with this pathway.\n");
		if (background != null)
			writer.write("# Pathway size in background: Intersection of genes in this pathway and the given background.\n");
		writer.write("# Pathway size: Number of genes in this pathway.\n");
		writer.write("# Molecules contributed enrichment: Names of query genes that overlaps with this pathway.\n");
		writer.write("Pathway name\tResource\tPathway Commons ID\tP-value\tQ-value\tHit size\t" + (background == null ? "" : "Pathway size in background\t") + "Pathway size\tMolecules contributed enrichment");
		for (String id : ids)
		{
			if (pvals[0].get(id) > 0.05) break;

			Set<String> all = new HashSet<>();
			if (hasGene) all.addAll(getGenes(id));
			if (hasChem) all.addAll(getChems(id));
			int allSize = all.size();
			Set<String> g = new HashSet<>(all);
			g.retainAll(molecules);
			int hitSize = g.size();
			int inBg = 0;
			if (background != null)
			{
				all.retainAll(background);
				inBg = all.size();
			}
			writer.write("\n" + getName(id) + "\t" + getResource(id) + "\t" + id + "\t" + pvals[0].get(id) + "\t" +
			qvals.get(id) + "\t" + hitSize + "\t" + (background == null ? "" : inBg + "\t") + allSize + "\t" + getMolsInString(g));
		}
		writer.close();
	}

	private String getMolsInString(Set<String> mols)
	{
		if (mols.isEmpty()) return "";
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for (String mol : mols)
		{
			if (first) first = false;
			else sb.append(", ");

			String name = ChEBI.get().getName(mol);
			sb.append(name == null ? mol : name);
		}
		return sb.toString();
	}

	/**
	 * Gets pathways sorted to their containment of the query genes. Does not control the density,
	 * but only the count, so it is likely that first pathways will be big ones.
	 */
	public TreeMap<String, Integer> getSortedPathways(Collection<String> genes)
	{
		final Map<String, Integer> map = new HashMap<>();

		for (String gene : genes)
		{
			for (String pathway : getPathways(gene))
			{
				if (map.containsKey(pathway)) map.put(pathway, map.get(pathway) + 1);
				else map.put(pathway, 1);
			}
		}

		TreeMap<String, Integer> sorted = new TreeMap<>((o1, o2) -> {
			return map.get(o2).compareTo(map.get(o1));
		});

		sorted.putAll(map);
		return sorted;
	}

	public static void main(String[] args) throws IOException
	{
//		String s = "PTEN, PIK3CA, ARID1A, PIK3R1, CTNNB1, TP53, KRAS, CTCF, FBXW7, LRP2, FGFR2, RYR1, TBL1XR1, MTOR, CACNA1A, PPP2R1A, PKN1, LYST, TRPM6, ERBB2, FN1, WDFY3, MYC, SPTB, DVL3, PRKCI, ECT2, ACTL6A, TBC1D31, IKBKB, PRKACA, DLG1, PTK2, THPO, DNM2, FOSL2, DSTYK, CCNE1, TNK2, EFNA1, PAK2, RASAL1, ARMC6, HGS, CDC37, TNFSF10, PPP1R1B, GRB2, PPP1CA";
//		List<String> select = getEnrichedPathways(Arrays.asList(s.split(", ")), null, 0.01);
//		for (String id : select)
//		{
//			System.out.println(id + "\t" + getName(id));
//		}

		prepareResource();
	}

	// Section: Preparing the pathways file from PC owl

	private static void prepareResource() throws IOException
	{
		Kronometre k = new Kronometre();
		Map<String, Set<String>> pathway2gene = new ConcurrentHashMap<>();
		Map<String, Set<String>> pathway2chem = new ConcurrentHashMap<>();
		Map<String, String> pathway2name = new ConcurrentHashMap<>();
		Map<String, String> pathway2resource = new ConcurrentHashMap<>();

		SimpleIOHandler h = new SimpleIOHandler();
		Model model = h.convertFromOWL(new FileInputStream("/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl"));
		Blacklist blacklist = new Blacklist("/home/babur/Documents/PC/blacklist.txt");

		System.out.println("Loaded BioPAX model");
		k.print();
		k.start();

		model.getObjects(Pathway.class).parallelStream().forEach(pathway ->
		{
			String id = pathway.getUri();
			String name = pathway.getDisplayName();

			if (name == null || name.isEmpty()) return;
			if (name.contains(" = ") || name.contains(" => ") || name.contains("[[") || name.contains("]]")) return;

			pathway2name.put(id, name);

			String resource = pathway.getDataSource().iterator().next().getName().iterator().next();
			pathway2resource.put(id, resource);

			Model m = excise(model, pathway);
			Set<String> syms = collectGeneSymbols(m);
			Set<String> chems = collectChemicals(m, blacklist);

			if (syms.size() + chems.size() < 3) return;

			if (!pathway2gene.containsKey(id)) pathway2gene.put(id, new HashSet<>());
			syms.forEach(pathway2gene.get(id)::add);

			if (!pathway2chem.containsKey(id)) pathway2chem.put(id, new HashSet<>());
			chems.forEach(pathway2chem.get(id)::add);
		});
		writeResources(pathway2gene, pathway2chem, pathway2name, pathway2resource);
		k.print();
	}

	private static void writeResources(Map<String, Set<String>> pathway2gene, Map<String, Set<String>> pathway2chem,
		Map<String, String> pathway2name, Map<String, String> pathway2resource) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("../repo/resource-files/" + FILE));

		for (String id : pathway2gene.keySet())
		{
			if (pathway2gene.get(id).isEmpty()) continue;
			if (pathway2name.get(id).contains("$")) continue;

			writer.write(id + "\t" + pathway2name.get(id) + "\t" + pathway2resource.get(id));

			for (String sym : pathway2gene.get(id))
			{
				writer.write("\t" + sym);
			}

			for (String chem : pathway2chem.get(id))
			{
				writer.write("\t" + chem);
			}

			writer.write("\n");
		}
		writer.close();
	}

	private static Model excise(Model model, Pathway pathway)
	{
		Completer c = new Completer(SimpleEditorMap.L3);

		Set<BioPAXElement> objects = c.complete(Collections.<BioPAXElement>singleton(pathway), model);

		Cloner cln = new Cloner(SimpleEditorMap.L3, BioPAXLevel.L3.getDefaultFactory());

		return cln.clone(model, objects);
	}

	private static Set<String> collectGeneSymbols(Model model)
	{
		return model.getObjects(EntityReference.class).stream()
			.map(EntityReference::getXref)
			.flatMap(Collection::stream)
			.filter(x -> x.getDb() != null && x.getDb().equalsIgnoreCase("HGNC Symbol"))
			.map(Xref::getId)
			.filter(Objects::nonNull)
			.map(id -> HGNC.get().getSymbol(id))
			.filter(Objects::nonNull)
			.collect(Collectors.toSet());
	}

	private static Set<String> collectChemicals(Model model, Blacklist blacklist)
	{
		return model.getObjects(EntityReference.class).stream()
			.filter(er -> !blacklist.getListed().contains(er.getUri()))
			.map(EntityReference::getXref)
			.flatMap(Collection::stream)
			.filter(x -> x.getDb() != null && x.getDb().equalsIgnoreCase("ChEBI"))
			.map(Xref::getId)
			.filter(Objects::nonNull)
			.filter(id -> id.startsWith("CHEBI:"))
			.collect(Collectors.toSet());
	}
}
