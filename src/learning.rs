struct Object {
    id: usize,
    parent_id: Option<usize>, // Assuming parent ID might not exist for some objects
    children: Vec<usize>,    // IDs of child objects
}

impl Object {
    fn new(id: usize, parent_id: Option<usize>) -> Self {
        Self {
            id,
            parent_id,
            children: Vec::new(),
        }
    }
}

struct ObjectManager {
    objects: Vec<Object>,
}

impl ObjectManager {
    fn build_hierarchy(&mut self) {
        // Temporary vector to store (parent_id, child_id) relationships
        let mut updates: Vec<(usize, usize)> = Vec::new();

        // First pass: Collect parent-child relationships
        for obj in &self.objects {
            if let Some(parent_id) = obj.parent_id {
                updates.push((parent_id, obj.id));
            }
        }

        // Second pass: Apply updates to parent objects
        for (parent_id, child_id) in updates {
            if let Some(parent_obj) = self.objects.get_mut(parent_id) {
                parent_obj.children.push(child_id);
            } else {
                eprintln!("Parent object with ID {} not found", parent_id);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use bio::utils::Interval;
    use super::*;

    #[test]
    fn test_object_mod() {
        let mut manager = ObjectManager {
            objects: vec![
                Object::new(0, None),            // Root object
                Object::new(1, Some(0)),        // Child of object 0
                Object::new(2, Some(0)),        // Child of object 0
                Object::new(3, Some(1)),        // Child of object 1
            ],
        };
    
        manager.build_hierarchy();
    
        for obj in &manager.objects {
            println!(
                "Object ID: {}, Parent ID: {:?}, Children: {:?}",
                obj.id, obj.parent_id, obj.children
            );
        }
    }

    #[test]
    fn test_mut_interval() {
        let mut interval = Interval::new(1..10).unwrap();
        assert_eq!(interval.start, 1);
        assert_eq!(interval.end, 10);
        interval.start = 5;
        interval.end = 15;
        assert_eq!(interval.start, 5);
        assert_eq!(interval.end, 15);
    }
}