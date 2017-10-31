package org.panda.resource;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Provides GO ontology graph.
 *
 * @author Ozgun Babur
 */
public class GO extends FileServer
{
	private Map<String, String> idToName;
	private Map<String, Set<String>> geneToGO;
	private Map<String, Set<String>> goToGene;
	private Map<String, Set<String>> isAMap;

	private static GO instance;

	public String getNameOfTerm(String goID)
	{
		return idToName.get(goID);
	}

	public Set<String> getGOIDs(String gene)
	{
		if (geneToGO.containsKey(gene)) return geneToGO.get(gene);
		return Collections.emptySet();
	}

	public Set<String> getAllGenes()
	{
		return geneToGO.keySet();
	}

	public Set<String> getGenes(Set<String> terms)
	{
		Set<String> genes = new HashSet<>();
		for (String term : terms)
		{
			genes.addAll(getGenes(term));
		}
		return genes;
	}

	public Set<String> getGenes(String term)
	{
		if (goToGene.containsKey(term)) return goToGene.get(term);
		return Collections.emptySet();
	}

	public Set<String> getParentTerms(String term)
	{
		if (isAMap.containsKey(term)) return isAMap.get(term);
		return Collections.emptySet();
	}

	public void printAssociatedTerms(Set<String> genes)
	{
		printAssociatedTerms(genes, Collections.emptySet(), Collections.emptySet());
	}

	public void printAssociatedTerms(Set<String> genes, Set<String> ignoreTerms, Set<String> ignoreGenesOfTerms)
	{
		Set<String> termSet = genes.stream().map(this::getGOIDs).flatMap(Collection::stream).collect(Collectors.toSet());
		termSet.removeAll(ignoreTerms);

		genes = new HashSet<>(genes);
		genes.removeAll(getGenes(ignoreGenesOfTerms));

		Map<String, Set<String>> map = new HashMap<>();
		for (String term : termSet)
		{
			Set<String> set = new HashSet<>(getGenes(term));
			set.retainAll(genes);
			assert !set.isEmpty();
			map.put(term, set);
		}
		termSet.stream().sorted((t1, t2) -> new Integer(map.get(t2).size()).compareTo(map.get(t1).size()))
			.forEach(term -> System.out.println(term + "\t" + getNameOfTerm(term) + "\t" + map.get(term)));
	}

	public void printAssociatedTerms(String gene, Set<String> ignore)
	{
		System.out.println("\ngene = " + gene);

		if (geneToGO.containsKey(gene))
		{
			for (String term : geneToGO.get(gene))
			{
				if (ignore == null || !ignore.contains(term))
				{
					String name = idToName.get(term);

					if (!name.startsWith("positive regulation ") && !name.startsWith("negative regulation "))
					{
						System.out.println(term + "\t" + name);
					}
				}
			}
		}
	}

	public Map<String, String> getTermsContaining(String query)
	{
		return idToName.keySet().stream().filter(id -> idToName.get(id).contains(query))
			.collect(Collectors.toMap(id -> id, id -> idToName.get(id)));
	}

	public Set<String> getGenesContainingKeywordInTermNames(String keyword)
	{
		return geneToGO.keySet().stream().filter(g -> geneToGO.get(g).stream().map(idToName::get)
			.anyMatch(name -> name.contains(keyword))).collect(Collectors.toSet());
	}

	public Set<String> getGenesContainingKeywordInTermNames(Set<String> keywords)
	{
		Set<String> genes = new HashSet<>();
		for (String keyword : keywords)
		{
			genes.addAll(geneToGO.keySet().stream().filter(g -> geneToGO.get(g).stream().map(idToName::get)
				.anyMatch(name -> name.contains(keyword))).collect(Collectors.toSet()));
		}
		return genes;
	}

	public static GO get()
	{
		if (instance == null) instance = new GO();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"GO-ontology.txt", "GO-gene-association.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://purl.obolibrary.org/obo/go/go-basic.obo",
			"http://geneontology.org/gene-associations/goa_human.gaf.gz"};
	}

	@Override
	public boolean load() throws IOException
	{
		idToName = new HashMap<>();
		isAMap = new HashMap<>();
		Scanner sc = new Scanner(new File(ResourceDirectory.get() + File.separator + getLocalFilenames()[0]));
		while (sc.hasNextLine())
		{
			String line = sc.nextLine();
			if (line.startsWith("id: "))
			{
				String term = line.substring(4);
				idToName.put(term, sc.nextLine().substring(6));

				while (sc.hasNextLine() && !line.startsWith("["))
				{
					line = sc.nextLine();
					if (line.startsWith("is_a: "))
					{
						String parent = line.substring(6);
						parent = parent.substring(0, parent.indexOf(" "));
						if (!isAMap.containsKey(term)) isAMap.put(term, new HashSet<>());
						isAMap.get(term).add(parent);
					}
				}
			}
		}

		geneToGO = new HashMap<>();
		goToGene = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[1]).filter(l -> !l.startsWith("!")).map(l -> l.split("\t")).forEach(t -> {
			String gene = t[2];
			String rel = t[3];
			String go = t[4];

			if (rel.contains("NOT")) return;

			if (!geneToGO.containsKey(gene)) geneToGO.put(gene, new HashSet<>());
			geneToGO.get(gene).add(go);
			if (!goToGene.containsKey(go)) goToGene.put(go, new HashSet<>());
			goToGene.get(go).add(gene);
		});

		for (String term : new HashSet<>(goToGene.keySet()))
		{
			Set<String> genes = goToGene.get(term);

			for (String parent : getParentTerms(term))
			{
				addParentsRecursive(parent, genes);
			}
		}

		return true;
	}

	private void addParentsRecursive(String parent, Set<String> genes)
	{
		if (!goToGene.containsKey(parent)) goToGene.put(parent, new HashSet<>());
		goToGene.get(parent).addAll(genes);
		for (String gene : genes)
		{
			geneToGO.get(gene).add(parent);
		}

		for (String grand : getParentTerms(parent))
		{
			addParentsRecursive(grand, genes);
		}
	}

	// Section: Temporary methods

	public static void main(String[] args) throws IOException
	{
//		List<String> genes = readPanCanGenes();
//		Set<String> ignore = readIgnoreList();

//		get().printAssociatedTerms(new HashSet<>(genes));

//		genes.forEach(g -> get().printAssociatedTerms(g, ignore));
//		get().printAssociatedTerms("ZNF732", Collections.emptySet());

		// Print terms containing a keyword
//		Map<String, String> termMap = get().getTermsContaining("damage");
//		for (String id : termMap.keySet())
//		{
//			System.out.println(id + "\t" + termMap.get(id));
//		}

//		checkInterestTermAssociationRate();

		get().printAssociatedTerms("MPHOSPH10", null);
	}

	private static List<String> readPanCanGenes() throws IOException
	{
		return readFirstColumn("/home/babur/Documents/PanCan/pancan.txt", 20);
	}
	private static Set<String> readIgnoreList() throws IOException
	{
		return new HashSet<>(readFirstColumn("/home/babur/Documents/PanCan/GO-term-ignore-list.txt", Integer.MAX_VALUE));
	}

	private static List<String> readFirstColumn(String file, int lineLimit) throws IOException
	{
		return Files.lines(Paths.get(file)).skip(1).limit(lineLimit)
			.map(l -> l.split("\t")[0]).collect(Collectors.toList());
	}

	private static void checkInterestTermAssociationRate() throws IOException
	{
		Set<String> interest = Files.lines(Paths.get("/home/babur/Documents/PanCan/GO-term-interest.txt"))
			.filter(l -> !l.startsWith("#")).collect(Collectors.toSet());

		System.out.println("Total genes = " + get().geneToGO.keySet().size());

		for (String gene : get().geneToGO.keySet())
		{
			for (String id : get().geneToGO.get(gene))
			{
				if (get().idToName.get(id) == null) System.err.println("Hell!");
			}
		}

		Set<String> match = new HashSet<>();
		for (String s : interest)
		{
			Set<String> genes = get().getGenesContainingKeywordInTermNames(s);
			int cnt = genes.size();
			System.out.println(cnt + "\t" + s);
			match.addAll(genes);
		}
		System.out.println();
		System.out.println("match.size() = " + match.size());
	}
}
