//! Preset registry and configuration types.
//!
//! All biological knowledge enters through YAML preset files (`presets/`).
//! The engine implements the math that configuration parameterizes.
//! No biological constants are hardcoded in Rust — the YAML files are the
//! shared contract between SeqChain, Needletail, and GenomeHub.
//!
//! Presets are organized by category (axis of variation):
//!   - `recognition` — how to find target sites (PAM patterns, motifs)
//!   - `tiling` — how to partition the genome into regions (organism-specific)
//!
//! The [`PresetRegistry`] scans all embedded YAMLs and indexes them by
//! `(category, id)`. The method catalog and API routes query the registry
//! to programmatically generate select options for GenomeHub.

use std::collections::HashMap;

use serde::Deserialize;
use serde_json::{json, Value};

// ---------------------------------------------------------------------------
// Anchor — named landmark on a gene
// ---------------------------------------------------------------------------

/// Named landmark on a gene for anchor-based annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Anchor {
    FivePrime,
    ThreePrime,
    Midpoint,
    None,
}

impl Anchor {
    pub fn from_str_anchor(s: &str) -> Self {
        match s {
            "five_prime" => Anchor::FivePrime,
            "three_prime" => Anchor::ThreePrime,
            "midpoint" => Anchor::Midpoint,
            _ => Anchor::None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Anchor::FivePrime => "five_prime",
            Anchor::ThreePrime => "three_prime",
            Anchor::Midpoint => "midpoint",
            Anchor::None => "",
        }
    }
}

// ---------------------------------------------------------------------------
// CRISPRPreset — loaded from YAML
// ---------------------------------------------------------------------------

/// Configuration for a CRISPR nuclease system.
#[derive(Debug, Clone)]
pub struct CRISPRPreset {
    pub name: String,
    pub pam: String,
    pub spacer_len: usize,
    pub pam_direction: String,
    pub description: String,
    pub mismatches: u8,
}

/// Raw YAML shape for CRISPR presets.
#[derive(Deserialize)]
struct CrisprYaml {
    name: String,
    pam: String,
    spacer_len: usize,
    pam_direction: String,
    description: String,
}

impl From<CrisprYaml> for CRISPRPreset {
    fn from(y: CrisprYaml) -> Self {
        CRISPRPreset {
            name: y.name,
            pam: y.pam,
            spacer_len: y.spacer_len,
            pam_direction: y.pam_direction,
            description: y.description.trim().to_string(),
            mismatches: 2,
        }
    }
}

impl CRISPRPreset {
    /// Look up a preset by name. Case-insensitive.
    pub fn by_name(name: &str) -> Option<Self> {
        let key = name.to_lowercase();
        let key = if key == "cpf1" { "cas12a" } else { &key };
        PresetRegistry::yaml("recognition", key)
            .and_then(|yaml| serde_yaml::from_str::<CrisprYaml>(yaml).ok())
            .map(CRISPRPreset::from)
    }

    /// List all available preset names.
    pub fn list() -> Vec<&'static str> {
        PresetRegistry::list("recognition")
    }
}

// ---------------------------------------------------------------------------
// FeatureConfig — loaded from YAML
// ---------------------------------------------------------------------------

/// A feature type defined by spatial relation to genes.
#[derive(Debug, Clone)]
pub struct FeatureDefinition {
    pub name: String,
    pub relation: String,
    pub max_distance: i64,
    pub priority: i32,
    pub anchor: Anchor,
}

/// Configuration for feature-type annotation.
#[derive(Debug, Clone)]
pub struct FeatureConfig {
    pub organism: String,
    pub features: Vec<FeatureDefinition>,
    pub default_feature: String,
}

/// Raw YAML shape for a single feature definition.
#[derive(Deserialize)]
struct FeatureDefYaml {
    relation: String,
    #[serde(default)]
    max_distance: i64,
    priority: i32,
    #[serde(default)]
    anchor: String,
}

/// Raw YAML shape for feature configs.
#[derive(Deserialize)]
struct FeatureYaml {
    organism: String,
    features: HashMap<String, FeatureDefYaml>,
    default_feature: String,
    #[serde(default)]
    priority: Vec<String>,
}

impl From<FeatureYaml> for FeatureConfig {
    fn from(y: FeatureYaml) -> Self {
        // Use the priority list to determine order; fall back to priority field.
        let mut defs: Vec<FeatureDefinition> = y
            .features
            .into_iter()
            .map(|(name, def)| FeatureDefinition {
                name,
                relation: def.relation,
                max_distance: def.max_distance,
                priority: def.priority,
                anchor: Anchor::from_str_anchor(&def.anchor),
            })
            .collect();

        if !y.priority.is_empty() {
            // Sort by position in the priority list.
            defs.sort_by_key(|d| {
                y.priority
                    .iter()
                    .position(|p| p == &d.name)
                    .unwrap_or(usize::MAX)
            });
        } else {
            defs.sort_by_key(|d| d.priority);
        }

        FeatureConfig {
            organism: y.organism,
            features: defs,
            default_feature: y.default_feature,
        }
    }
}

impl FeatureConfig {
    /// Look up a feature config by name. Case-insensitive.
    pub fn by_name(name: &str) -> Option<Self> {
        let key = name.to_lowercase();
        let key = if key == "saccer3_features" { "saccer3" } else { &key };
        PresetRegistry::yaml("tiling", key)
            .and_then(|yaml| serde_yaml::from_str::<FeatureYaml>(yaml).ok())
            .map(FeatureConfig::from)
    }

    /// List all available config names.
    pub fn list() -> Vec<&'static str> {
        PresetRegistry::list("tiling")
    }
}

// ---------------------------------------------------------------------------
// PresetRegistry — unified index over all preset categories
// ---------------------------------------------------------------------------

/// A single embedded preset: its category, id, and raw YAML source.
struct EmbeddedPreset {
    category: &'static str,
    id: &'static str,
    yaml: &'static str,
}

/// All embedded presets. Adding a new YAML = one line here.
static ALL_PRESETS: &[EmbeddedPreset] = &[
    // recognition
    EmbeddedPreset { category: "recognition", id: "spcas9", yaml: include_str!("../../../../presets/recognition/spcas9.yaml") },
    EmbeddedPreset { category: "recognition", id: "cas12a", yaml: include_str!("../../../../presets/recognition/cas12a.yaml") },
    EmbeddedPreset { category: "recognition", id: "cas12m", yaml: include_str!("../../../../presets/recognition/cas12m.yaml") },
    // tiling
    EmbeddedPreset { category: "tiling", id: "saccer3", yaml: include_str!("../../../../presets/tiling/saccer3_features.yaml") },
];

/// Generic YAML value for extracting label/description from any preset.
#[derive(Deserialize)]
struct PresetMeta {
    #[serde(default)]
    name: Option<String>,
    #[serde(default)]
    organism: Option<String>,
    #[serde(default)]
    description: Option<String>,
}

impl PresetMeta {
    /// Best human-readable label for this preset.
    fn label(&self) -> String {
        self.name
            .as_deref()
            .or(self.organism.as_deref())
            .unwrap_or("?")
            .to_string()
    }

    fn description(&self) -> Option<String> {
        self.description.as_ref().map(|d| d.trim().to_string())
    }
}

/// Registry over all preset categories. Provides generic listing and
/// select-option generation so the method catalog and API routes never
/// need to know what categories exist.
pub struct PresetRegistry;

impl PresetRegistry {
    /// List all known categories.
    pub fn categories() -> Vec<&'static str> {
        let mut cats: Vec<&str> = ALL_PRESETS.iter().map(|p| p.category).collect();
        cats.dedup();
        cats
    }

    /// List preset IDs in a category.
    pub fn list(category: &str) -> Vec<&'static str> {
        ALL_PRESETS
            .iter()
            .filter(|p| p.category == category)
            .map(|p| p.id)
            .collect()
    }

    /// Get raw YAML for a preset.
    pub fn yaml(category: &str, id: &str) -> Option<&'static str> {
        ALL_PRESETS
            .iter()
            .find(|p| p.category == category && p.id == id)
            .map(|p| p.yaml)
    }

    /// Build select options for a category — ready to embed in the method catalog.
    ///
    /// Conforms to the GenomeHub engine interface contract:
    /// ```json
    /// {
    ///   "value": "id",
    ///   "label": "Human Name",
    ///   "description": "Prose description",
    ///   "parameters": { ... domain-specific fields ... }
    /// }
    /// ```
    ///
    /// GenomeHub renders `label` + `description` in the dropdown and displays
    /// `parameters` as read-only key-value pairs when the option is selected.
    pub fn select_options(category: &str) -> Vec<Value> {
        // Fields that are part of the select option envelope, not domain parameters.
        const ENVELOPE: &[&str] = &["name", "type", "description"];

        ALL_PRESETS
            .iter()
            .filter(|p| p.category == category)
            .filter_map(|p| {
                let val: Value = serde_yaml::from_str(p.yaml).ok()?;
                let meta: PresetMeta = serde_yaml::from_str(p.yaml).ok()?;
                let obj = val.as_object()?;

                // Domain-specific fields go into "parameters"
                let params: serde_json::Map<String, Value> = obj
                    .iter()
                    .filter(|(k, _)| !ENVELOPE.contains(&k.as_str()))
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect();

                Some(json!({
                    "value": p.id,
                    "label": meta.label(),
                    "description": meta.description().unwrap_or_default(),
                    "parameters": params,
                }))
            })
            .collect()
    }

    /// Full detail for all presets in a category (for `/api/presets/:category`).
    pub fn detail(category: &str) -> Vec<Value> {
        ALL_PRESETS
            .iter()
            .filter(|p| p.category == category)
            .filter_map(|p| {
                let mut val: Value = serde_yaml::from_str(p.yaml).ok()?;
                if let Some(obj) = val.as_object_mut() {
                    obj.insert("id".into(), json!(p.id));
                }
                Some(val)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_anchor_roundtrip() {
        for anchor in [
            Anchor::FivePrime,
            Anchor::ThreePrime,
            Anchor::Midpoint,
            Anchor::None,
        ] {
            assert_eq!(Anchor::from_str_anchor(anchor.as_str()), anchor);
        }
    }

    #[test]
    fn test_spcas9_preset() {
        let p = CRISPRPreset::by_name("spcas9").unwrap();
        assert_eq!(p.pam, "NGG");
        assert_eq!(p.spacer_len, 20);
        assert_eq!(p.pam_direction, "downstream");
    }

    #[test]
    fn test_cas12a_preset() {
        let p = CRISPRPreset::by_name("cas12a").unwrap();
        assert_eq!(p.pam, "TTTN");
        assert_eq!(p.spacer_len, 24);
        assert_eq!(p.pam_direction, "upstream");
    }

    #[test]
    fn test_cas12m_preset() {
        let p = CRISPRPreset::by_name("cas12m").unwrap();
        assert_eq!(p.pam, "TTN");
        assert_eq!(p.spacer_len, 20);
        assert_eq!(p.pam_direction, "downstream");
    }

    #[test]
    fn test_cpf1_alias() {
        let p = CRISPRPreset::by_name("cpf1").unwrap();
        assert_eq!(p.name, "Cas12a");
    }

    #[test]
    fn test_list_presets() {
        let names = CRISPRPreset::list();
        assert!(names.contains(&"spcas9"));
        assert!(names.contains(&"cas12a"));
        assert!(names.contains(&"cas12m"));
    }

    #[test]
    fn test_saccer3_config() {
        let c = FeatureConfig::by_name("saccer3").unwrap();
        assert_eq!(c.features.len(), 3);
        assert_eq!(c.features[0].name, "gene_body");
        assert_eq!(c.features[1].name, "promoter");
        assert_eq!(c.features[1].max_distance, 500);
        assert_eq!(c.default_feature, "intergenic");
    }

    #[test]
    fn test_saccer3_features_alias() {
        assert!(FeatureConfig::by_name("saccer3_features").is_some());
    }

    #[test]
    fn test_registry_categories() {
        let cats = PresetRegistry::categories();
        assert!(cats.contains(&"recognition"));
        assert!(cats.contains(&"tiling"));
    }

    #[test]
    fn test_registry_select_options() {
        let opts = PresetRegistry::select_options("recognition");
        assert_eq!(opts.len(), 3);
        assert_eq!(opts[0]["value"], "spcas9");
        assert!(opts[0]["label"].as_str().unwrap().contains("Cas9"));

        let opts = PresetRegistry::select_options("tiling");
        assert_eq!(opts.len(), 1);
        assert_eq!(opts[0]["value"], "saccer3");
    }

    #[test]
    fn test_registry_detail() {
        let detail = PresetRegistry::detail("recognition");
        assert_eq!(detail.len(), 3);
        // Each entry should have an injected "id" field
        assert_eq!(detail[0]["id"], "spcas9");
        assert!(detail[0]["pam"].is_string());
    }
}
